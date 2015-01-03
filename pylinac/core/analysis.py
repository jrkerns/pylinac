from abc import ABCMeta, abstractmethod


class AnalysisModule(metaclass=ABCMeta):
    """An abstract class for distinct analysis modules (Starshot, VMAT, etc). Its purpose is to define basic method
    names to be consistent across analysis modules. E.g. the starshot and VMAT modules will both have an "analyze" and
    "load_image" method.
    """
    @abstractmethod
    def load_demo_image(self):
        pass

    @abstractmethod
    def analyze(self):
        pass

    @abstractmethod
    def run_demo(self):
        pass

    @abstractmethod
    def load_image(self, filepath):
        pass

    @abstractmethod
    def load_image_UI(self):
        pass