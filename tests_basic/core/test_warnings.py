import warnings
from unittest import TestCase

from pylinac.core.warnings import WarningCollectorMixin, capture_warnings


class TestWarnings(TestCase):
    def test_warnings_are_captured(self):
        @capture_warnings
        class MyClass(WarningCollectorMixin):
            def my_method(self):
                warnings.warn("This is a warning", UserWarning)
                warnings.warn("This is another warning", DeprecationWarning)

        my_instance = MyClass()
        my_instance.my_method()

        self.assertEqual(len(my_instance.get_captured_warnings()), 2)
        first_warning = my_instance.get_captured_warnings()[0]
        self.assertEqual(first_warning["message"], "This is a warning")
        self.assertEqual(first_warning["category"], "UserWarning")
        second_warning = my_instance.get_captured_warnings()[1]
        self.assertEqual(second_warning["message"], "This is another warning")
        self.assertEqual(second_warning["category"], "DeprecationWarning")
