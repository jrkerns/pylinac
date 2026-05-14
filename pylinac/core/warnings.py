import sys
import types
import warnings
from functools import wraps
from threading import Lock
from typing import Any, TypeVar

T = TypeVar("T", bound=type[Any])  # Class type


class WarningCollectorMixin:
    """A mixin class that captures warnings emitted by its methods."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._captured_warnings: list[dict] = []
        self._warnings_lock: Lock = Lock()
        self._in_warning_capture: bool = False

    def _add_warnings(self, warnings_list: list[dict]) -> None:
        with self._warnings_lock:
            self._captured_warnings.extend(warnings_list)

    def get_captured_warnings(self) -> list[dict]:
        """Retrieve the list of captured warnings, deduplicated."""
        with self._warnings_lock:
            seen: set[tuple] = set()
            unique: list[dict] = []
            for w in self._captured_warnings:
                key = tuple(sorted(w.items()))
                if key not in seen:
                    seen.add(key)
                    unique.append(w)
            return unique

    def clear_captured_warnings(self) -> None:
        """Clear the list of captured warnings."""
        with self._warnings_lock:
            self._captured_warnings.clear()


def capture_warnings_method_wrapper(method: callable) -> callable:
    """Decorator to capture warnings emitted by the decorated method.

    All warning categories are captured and stored. They are also
    re-emitted to stderr so they remain visible in the console.

    A nesting guard (``_in_warning_capture``) ensures that when wrapped
    methods call other wrapped methods, only the outermost context
    captures warnings — avoiding duplicate entries.
    """

    @wraps(method)
    def wrapper(self, *args, **kwargs):
        if getattr(self, "_in_warning_capture", False):
            return method(self, *args, **kwargs)

        self._in_warning_capture = True
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = method(self, *args, **kwargs)
                captured = []
                for warning in w:
                    warning_info = {
                        "message": str(warning.message),
                        "category": warning.category.__name__,
                        "filename": warning.filename,
                        "lineno": warning.lineno,
                        "line": warning.line,
                    }
                    captured.append(warning_info)
            self._add_warnings(captured)
            for warning in w:
                warnings.showwarning(
                    message=warning.message,
                    category=warning.category,
                    filename=warning.filename,
                    lineno=warning.lineno,
                    file=sys.stderr,
                    line=warning.line,
                )
            return result
        finally:
            self._in_warning_capture = False

    return wrapper


def capture_warnings(cls: type[T]) -> type[T]:
    """Class decorator to apply warning capture to all public methods.

    Methods defined directly on ``cls`` are wrapped first.  Then the MRO
    is walked so that inherited public methods (e.g. ``analyze`` on a
    base class) are also wrapped on the decorated subclass.  Private and
    dunder methods from parent classes are skipped.
    """
    for attr_name, attr in cls.__dict__.items():
        if isinstance(attr, (types.FunctionType, types.MethodType)):
            setattr(cls, attr_name, capture_warnings_method_wrapper(attr))

    for parent in cls.__mro__[1:]:
        if parent is object:
            continue
        for attr_name, attr in parent.__dict__.items():
            if attr_name.startswith("_"):
                continue
            if attr_name in cls.__dict__:
                continue
            if isinstance(attr, (types.FunctionType, types.MethodType)):
                setattr(cls, attr_name, capture_warnings_method_wrapper(attr))
    return cls
