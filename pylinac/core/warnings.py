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

    def _add_warnings(self, warnings_list: list[dict]) -> None:
        with self._warnings_lock:
            self._captured_warnings.extend(warnings_list)

    def get_captured_warnings(self) -> list[dict]:
        """Retrieve the list of captured warnings."""
        with self._warnings_lock:
            return list(self._captured_warnings)

    def clear_captured_warnings(self) -> None:
        """Clear the list of captured warnings."""
        with self._warnings_lock:
            self._captured_warnings.clear()


def capture_warnings_method_wrapper(method: callable) -> callable:
    """Decorator to capture warnings emitted by the decorated method."""

    @wraps(method)
    def wrapper(self, *args, **kwargs):
        with warnings.catch_warnings(record=True) as w:
            # Ensure all warnings are captured
            warnings.simplefilter("always")
            # Execute the original method
            result = method(self, *args, **kwargs)
            # Process captured warnings
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
            # Add captured warnings to the mixin's list
        self._add_warnings(captured)
        # Re-emit the warning outside the capture context
        # so it goes to the console as well
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

    return wrapper


def capture_warnings(cls: type[T]) -> type[T]:
    """
    Class decorator to automatically apply the capture_warnings decorator
    to all methods of the class that are not private or special methods.
    """
    for attr_name, attr in cls.__dict__.items():
        # Only decorate callable attributes (methods)
        if isinstance(attr, (types.FunctionType, types.MethodType)):
            setattr(cls, attr_name, capture_warnings_method_wrapper(attr))
    return cls
