

class AnalysisModule(object):
    """An abstract class for distinct analysis modules (Starshot, VMAT, etc). Its purpose is to define basic method
    names to be consistent across analysis modules. E.g. the starshot and VMAT modules will both have an "analyze" and
    "load_image" method.
    """

    def load_demo_image(self):
        """To be overloaded by each specific tool. Loads a demo image for the given class."""
        raise NotImplementedError("Loading the demo image(s) for this module has not been implemented yet.")

    def analyze(self):
        """To be overloaded by subclass. Main analysis method for module"""
        raise NotImplementedError("Analysis has not been implemented for this module.")

    def run_demo(self):
        """Demo of module's abilities."""
        raise NotImplementedError("The demo for this module has not been built yet.")

    def load_image(self, filepath):
        raise NotImplementedError

    def load_image_UI(self):
        raise NotImplementedError