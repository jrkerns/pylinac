import warnings
from unittest import TestCase

from pylinac.core.warnings import WarningCollectorMixin, capture_warnings


class TestWarnings(TestCase):
    def test_all_warning_categories_are_captured(self):
        @capture_warnings
        class MyClass(WarningCollectorMixin):
            def my_method(self):
                warnings.warn("This is a warning", UserWarning)
                warnings.warn("This is another warning", DeprecationWarning)

        my_instance = MyClass()
        my_instance.my_method()

        captured = my_instance.get_captured_warnings()
        self.assertEqual(len(captured), 2)
        self.assertEqual(captured[0]["message"], "This is a warning")
        self.assertEqual(captured[0]["category"], "UserWarning")
        self.assertEqual(captured[1]["message"], "This is another warning")
        self.assertEqual(captured[1]["category"], "DeprecationWarning")

    def test_inherited_method_warnings_are_captured(self):
        """Warnings from methods defined on a base class should be captured
        when the subclass is decorated with @capture_warnings."""

        class Base(WarningCollectorMixin):
            def analyze(self):
                self._do_work()

            def _do_work(self):
                warnings.warn(
                    "Could not determine phantom roll. Setting roll to 0.",
                    UserWarning,
                )

        @capture_warnings
        class Child(Base):
            pass

        instance = Child()
        instance.analyze()

        captured = instance.get_captured_warnings()
        self.assertEqual(len(captured), 1)
        self.assertEqual(
            captured[0]["message"],
            "Could not determine phantom roll. Setting roll to 0.",
        )

    def test_nested_wrapped_methods_do_not_duplicate(self):
        """When a wrapped method calls another wrapped method, warnings
        should only be captured once by the outermost context."""

        class Base(WarningCollectorMixin):
            def analyze(self):
                self.inner()

            def inner(self):
                warnings.warn("inner warning", UserWarning)

        @capture_warnings
        class Child(Base):
            pass

        instance = Child()
        instance.analyze()

        captured = instance.get_captured_warnings()
        self.assertEqual(len(captured), 1)
        self.assertEqual(captured[0]["message"], "inner warning")

    def test_duplicate_warnings_are_deduplicated(self):
        """Identical warnings emitted multiple times should appear only once
        in the output of get_captured_warnings()."""

        @capture_warnings
        class MyClass(WarningCollectorMixin):
            def run(self):
                for _ in range(5):
                    warnings.warn("repeated warning", UserWarning)
                warnings.warn("unique warning", RuntimeWarning)

        instance = MyClass()
        instance.run()

        captured = instance.get_captured_warnings()
        messages = [w["message"] for w in captured]
        self.assertEqual(messages.count("repeated warning"), 1)
        self.assertEqual(messages.count("unique warning"), 1)
        self.assertEqual(len(captured), 2)

    def test_caller_can_filter_reemitted_warning(self):
        @capture_warnings
        class MyClass(WarningCollectorMixin):
            def run(self):
                warnings.warn("bbox_area is deprecated", UserWarning)
                warnings.warn("another warning", UserWarning)

        instance = MyClass()
        with warnings.catch_warnings(record=True) as emitted:
            warnings.simplefilter("always")
            warnings.filterwarnings(
                "ignore", message=".*bbox_area.*", category=UserWarning
            )
            instance.run()

        self.assertEqual([str(w.message) for w in emitted], ["another warning"])
        self.assertEqual(
            [w["message"] for w in instance.get_captured_warnings()],
            ["another warning"],
        )

    def test_caller_can_ignore_warning_before_capture(self):
        @capture_warnings
        class MyClass(WarningCollectorMixin):
            def run(self):
                warnings.warn("bbox_area is deprecated", UserWarning)

        instance = MyClass()
        with warnings.catch_warnings(record=True) as emitted:
            warnings.filterwarnings(
                "ignore", message=".*bbox_area.*", category=UserWarning
            )
            instance.run()

        self.assertEqual(emitted, [])
        self.assertEqual(instance.get_captured_warnings(), [])

    def test_repeat_suppression_across_method_calls(self):
        @capture_warnings
        class MyClass(WarningCollectorMixin):
            def run(self):
                warnings.warn("repeated across calls", UserWarning)

        for action in ("once", "default", "module"):
            with self.subTest(action=action):
                instance = MyClass()
                with warnings.catch_warnings(record=True) as emitted:
                    warnings.simplefilter(action)
                    instance.run()
                    instance.run()

                self.assertEqual(len(emitted), 1)
                self.assertEqual(len(instance.get_captured_warnings()), 1)

    def test_always_filter_emits_on_every_call(self):
        @capture_warnings
        class MyClass(WarningCollectorMixin):
            def run(self):
                warnings.warn("repeated across calls", UserWarning)

        instance = MyClass()
        with warnings.catch_warnings(record=True) as emitted:
            warnings.simplefilter("always")
            instance.run()
            instance.run()

        self.assertEqual(len(emitted), 2)
        self.assertEqual(len(instance.get_captured_warnings()), 1)
