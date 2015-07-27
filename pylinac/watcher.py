"""The watcher file is a script meant to be run as an ongoing process to watch a given directory and analyzing files
that may be moved there for certain keywords. Automatic processing will be started if the file contains the keywords."""

import time
import logging
import os.path as osp
import sys
import abc

try:
    from watchdog.observers import Observer
    from watchdog.events import FileSystemEventHandler
except ImportError:
    raise ImportError("Watchdog must be installed to perform file watching.")

from pylinac.starshot import Starshot
from pylinac.picketfence import PicketFence
from pylinac.cbct import CBCT
from pylinac.log_analyzer import MachineLog
from pylinac.vmat import VMAT

logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')


class AnalyzeMixin:
    """Mixin for processing files caught by the pylinac watcher.

    Attributes
    ----------
    obj : class
        The class that analyzes the file; e.g. Starshot, PicketFence, etc.
    keywords : iterable
        Holds the keywords that are looked for in the file name. If the filename has that keyword, analysis will be attempted on that file.
    args : dict
        Dictionary that holds the tolerance settings of analysis.
    """
    obj = 'class that analyzes'
    keywords = ('',)
    save_image_method = 'save_analyzed_image'
    save_text_method = 'return_results'
    args = {}

    def __init__(self, path):
        self.basepath = osp.splitext(path)[0]
        self.instance = self.obj(path)
        self.analyze()
        self.save_image()
        self.save_text()

    @property
    def img_filename(self):
        """The name of the file for the analyzed image."""
        return self.basepath + '.png'

    @property
    def txt_filename(self):
        """The name of the file for the text results."""
        return self.basepath + '.txt'

    def save_image(self):
        """Save the analyzed image to file."""
        method = getattr(self.instance, self.save_image_method)
        method(self.img_filename)

    def save_text(self):
        """Save the analysis results to a text file."""
        method = getattr(self.instance, self.save_text_method)
        with open(self.txt_filename, 'w') as txtfile:
            txtfile.write(method())

    @abc.abstractproperty
    def analyze(self):
        pass


class AnalyzeStar(AnalyzeMixin):
    """Analysis class for starshots."""
    obj = Starshot
    keywords = ('star',)
    args = {'tolerance': 1, 'radius': 0.8}

    def analyze(self):
        self.instance.analyze(tolerance=self.args['tolerance'], radius=self.args['radius'])


class AnalyzePF(AnalyzeMixin):
    """Analysis class for picket fences."""
    obj = PicketFence
    keywords = ('pf', 'picket')
    args = {'tolerance': 0.5, 'action_tolerance': 0.3}

    def analyze(self):
        self.instance.analyze(tolerance=self.args['tolerance'], action_tolerance=self.args['action_tolerance'])


class AnalyzeCBCT(AnalyzeMixin):
    """Analysis class for CBCTs."""
    obj = CBCT.from_zip_file
    keywords = ('cbct', 'ct')
    args = {'hu tolerance': 40, 'scaling tolerance': 1}

    def analyze(self):
        self.instance.analyze(hu_tolerance=self.args['hu tolerance'], scaling_tolerance=self.args['scaling tolerance'])


class AnalyzeVMAT(AnalyzeMixin):
    """Analysis class for VMATs."""
    obj = VMAT.from_zip
    keywords = ('vmat', 'drgs', 'drmlc')
    args = {'tolerance': 1.5}

    def analyze(self):
        self.instance.analyze(tolerance=self.args['tolerance'])


class AnalyzeLog(AnalyzeMixin):
    """Analysis class for dynalogs or trajectory logs."""
    obj = MachineLog
    keywords = ('',)
    args = {'resolution': 0.1, 'distTA': 1, 'doseTA': 1, 'threshold': 10}
    save_image_method = 'save_summary'

    def analyze(self):
        self.instance.fluence.gamma.calc_map(doseTA=self.args['doseTA'], distTA=self.args['distTA'],
                                             resolution=self.args['resolution'], threshold=self.args['threshold'])

    def save_text(self):
        with open(self.txt_filename, 'w') as txtfile:
            txtfile.write(self.instance.report_basic_parameters(False))


def analysis_should_be_done(path):
    """Return boolean of whether the file should be analysed, based on if the filename
    has a keyword."""
    path = osp.basename(path).lower()
    for analysis_class in (AnalyzeStar, AnalyzeCBCT, AnalyzeVMAT, AnalyzeLog, AnalyzePF):
        for keyword in analysis_class.keywords:
            if (keyword in path) and not any(item in path for item in ('.png', '.txt')):
                # more specific filtering of data by type
                if analysis_class in (AnalyzeCBCT, AnalyzeVMAT):
                    if path.endswith('.zip'):
                        time.sleep(2)
                        return True, analysis_class
                elif analysis_class == AnalyzeLog:
                    if path.endswith('.dlg') or path.endswith('.bin'):
                        time.sleep(1)
                        return True, analysis_class
                else:
                    time.sleep(1)
                    return True, analysis_class
    return False, None


class FileAnalyzerEvent(FileSystemEventHandler):
    """Handler for file events."""

    def on_created(self, event):
        full_file_path = osp.abspath(event.src_path)
        # determine if file has keyword in it
        do_analysis, analysis_class = analysis_should_be_done(full_file_path)
        # process it if so
        if do_analysis:
            logging.info(full_file_path + " file found and will be analyzed...")
            analysis_class(full_file_path)
            logging.info(full_file_path + " was analyzed and now has an associated .txt and .png file")
        else:
            logging.info(full_file_path + " was added but was not deemed a file to be analyzed.")


if __name__ == "__main__":
    event_handler = FileAnalyzerEvent()
    observer = Observer()
    path = sys.argv[1] if len(sys.argv) > 1 else '.'
    observer.schedule(event_handler, path, recursive=True)
    logging.info("Pylinac now watching at " + osp.abspath(path))
    observer.start()
    while True:
        time.sleep(1)
    observer.join()


