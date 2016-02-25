"""The watcher file is a script meant to be run as an ongoing process to watch a given directory and analyzing files
that may be moved there for certain keywords. Automatic processing will be started if the file contains the keywords."""
import datetime
import logging
import os.path as osp
import time

try:
    from watchdog.observers import Observer
    from watchdog.events import FileSystemEventHandler
except ImportError:
    raise ImportError("Watchdog must be installed to perform file watching. Run ``pip install watchdog`` and try again.")
try:
    import yaml
except ImportError:
    raise ImportError("PyYaml must be installed to perform file watching. Run ``pip install pyyaml`` and try again.")

from pylinac import CBCT, VMAT, Starshot, PicketFence, MachineLog, WinstonLutz, LeedsTOR, PipsProQC3

logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')


class AnalyzeMixin:
    """Mixin for processing files caught by the pylinac watcher.

    Attributes
    ----------
    obj : class
        The class that analyzes the file; e.g. Starshot, PicketFence, etc.
    config_name : str
        The string that references the class in the YAML config file.
    save_image_method : str
        String-ified method name for saving the image to file.
    save_test_method : str
        String-ified method name for saving the text results to file.
    expecting_zip : bool
        Whether the class expects to find a ZIP archive, or normal file.
    """
    obj = object
    config_name = ''
    save_image_method = 'save_analyzed_image'
    save_text_method = 'return_results'
    expecting_zip = False

    def __init__(self, path, config):
        """
        Parameters
        ----------
        path : str
            The path to the file to be analyzed.
        config :
            The configuration settings of analysis. See `~pylinac.watcher.load_config`.
        """
        self.full_path = path
        self.local_path = osp.basename(path)
        self.base_name = osp.splitext(self.local_path)[0]
        self.config = config

    def process(self):
        """Process the file; includes analysis, saving results to file, and sending emails."""
        logging.info(self.local_path + " file found and will be analyzed...")
        self.instance = self.obj(self.full_path)
        self.analyze()
        self.save_image()
        self.save_text()
        if self.config['email']['enable-all']:
            self.send_email()
        elif self.config['email']['enable-failure'] and self.should_send_failure_email():
            self.send_email()
        logging.info(self.local_path + " was analyzed and now has an associated .txt and .png file")

    @property
    def img_filename(self):
        """The name of the file for the analyzed image."""
        return self.base_name + self.config['general']['file-suffix'] + '.png'

    @property
    def txt_filename(self):
        """The name of the file for the text results."""
        return self.base_name + self.config['general']['file-suffix'] + '.txt'

    @property
    def keywords(self):
        """The keywords that signal a file is of a certain analysis type."""
        return self.config[self.config_name]['keywords']

    def keyword_in_here(self):
        """Determine whether a keyword exists in the filename."""
        if not self.expecting_zip:
            return any(keyword in self.local_path.lower() for keyword in self.keywords)
        else:
            return any(keyword in self.local_path.lower() for keyword in self.keywords) and self.local_path.endswith('.zip')

    @property
    def failure_settings(self):
        """The YAML failure settings."""
        return self.config[self.config_name]['failure']

    @property
    def analysis_settings(self):
        """The YAML analysis settings."""
        return self.config[self.config_name]['analysis']

    def send_email(self):
        """Send an email with the analysis results."""
        try:
            import yagmail
        except ImportError:
            raise ImportError("Yagmail must be installed to perform file watching. Run ``pip install yagmail`` and try again.")

        # compose message
        current_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        if self.config['email']['enable-all']:
            statement = 'The pylinac watcher analyzed the file "{}" at {}. '
            statement += 'The analysis results are at "{}" and have also been attached here.'
        elif self.config['email']['enable-failure']:
            statement = 'The pylinac watcher analyzed the file "{}" at {} and '
            statement += 'found something that failed your configuration settings.'
            statement += 'The analysis results are at "{}" and have also been attached here.'
        statement = statement.format(self.local_path, current_time, osp.dirname(self.full_path))
        # send the email
        contents = [statement, self.img_filename, self.txt_filename]
        yagserver = yagmail.SMTP(self.config['email']['sender'], self.config['email']['sender-password'])
        recipients = [email for email in self.config['email']['recipients']]
        yagserver.send(to=recipients,
                       subject=self.config['email']['subject'],
                       contents=contents)
        logging.info("An email was sent to the recipients with the results")

    def save_image(self):
        """Save the analyzed image to file."""
        method = getattr(self.instance, self.save_image_method)
        method(self.img_filename)

    def save_text(self):
        """Save the analysis results to a text file."""
        method = getattr(self.instance, self.save_text_method)
        with open(self.txt_filename, 'w') as txtfile:
            txtfile.write(method())

    def should_send_failure_email(self):
        """Check whether analysis results were poor and an email should be triggered."""
        return not self.instance.passed

    def analyze(self):
        """Analyze the file."""
        self.instance.analyze(**self.analysis_settings)


class AnalyzeWL(AnalyzeMixin):
    """Analysis runner for Winston-Lutz images."""
    obj = WinstonLutz.from_zip
    save_text_method = 'results'
    save_image_method = 'save_summary'
    config_name = 'winston-lutz'
    expecting_zip = True

    def analyze(self):
        """Winson-Lutz doesn't explicitly analyze files."""
        pass

    def should_send_failure_email(self):
        """Failure of WL is based on 3 criteria."""
        send = False
        for key, val in self.failure_settings:
            if key == 'gantry-iso-size':
                if self.instance.gantry_iso_size > val:
                    send = True
            if key == 'mean-cax-bb-distance':
                if self.instance.cax2bb_distance() > val:
                    send = True
            if key == 'max-cax-bb-distance':
                if self.instance.cax2bb_distance('max') > val:
                    send = True
        return send


class AnalyzeLeeds(AnalyzeMixin):
    """Analysis runner for Leeds TOR phantom."""
    obj = LeedsTOR
    config_name = 'leeds'


class AnalyzePipsPro(AnalyzeMixin):
    """Analysis runner for PipsPro QC-3."""
    obj = PipsProQC3
    config_name = 'pipspro'


class AnalyzeStar(AnalyzeMixin):
    """Analysis runner for starshots."""
    obj = Starshot
    config_name = 'starshot'


class AnalyzePF(AnalyzeMixin):
    """Analysis runner for picket fences."""
    obj = PicketFence
    config_name = 'picketfence'


class AnalyzeCBCT(AnalyzeMixin):
    """Analysis runner for CBCTs."""
    obj = CBCT.from_zip
    config_name = 'cbct'
    expecting_zip = True

    def should_send_failure_email(self):
        """Failure of CBCT depends on individual module performance."""
        send = False
        for key, val in self.failure_settings:
            if key == 'hu-passed' and not self.instance.hu.overall_passed:
                send = True
            if key == 'uniformity-passed' and not self.instance.uniformity.overall_passed:
                send = True
            if key == 'geometry-passed' and not self.instance.geometry.overall_passed:
                send = True
            if key == 'thickness-passed' and not self.instance.thickness.passed:
                send = True
        return send


class AnalyzeVMAT(AnalyzeMixin):
    """Analysis runner for VMATs."""
    obj = VMAT.from_zip
    config_name = 'vmat'
    expecting_zip = True


class AnalyzeLog(AnalyzeMixin):
    """Analysis runner for dynalogs or trajectory logs."""
    obj = MachineLog
    save_image_method = 'save_summary'

    def analyze(self):
        """Log analysis is done via calculating gamma."""
        self.instance.fluence.gamma.calc_map(**self.analysis_settings)

    def keyword_in_here(self):
        """Log keywords just check file extensions."""
        return self.local_path.endswith('.dlg') or self.local_path.endswith('.bin')

    def save_text(self):
        """Special text save method."""
        with open(self.txt_filename, 'w') as txtfile:
            txtfile.write(self.instance.report_basic_parameters(False))

    def should_send_failure_email(self):
        """Failure is based on several varying criteria."""
        send = False
        for key, val in self.failure_settings:
            if key == 'gamma':
                gamma_below_threshold = self.instance.fluence.gamma.pass_prcnt < val
                if gamma_below_threshold:
                    send = True
            elif key == 'avg-rms':
                rms_above_threshold = self.instance.mlc.get_RMS_avg() > val
                if rms_above_threshold:
                    send = True
            elif key == 'max-rms':
                rms_above_threshold = self.instance.mlc.get_RMS_max() > val
                if rms_above_threshold:
                    send = True
        return send


def analysis_should_be_done(path, config):
    """Return boolean of whether the file should be analysed, based on if the filename has a keyword."""
    # return if not an analysis-worthy file
    if any(item in path for item in config['general']['avoid-keywords']):
        return False, None
    else:
        analysis_classes = (AnalyzeStar, AnalyzeCBCT, AnalyzeVMAT, AnalyzePF, AnalyzeWL, AnalyzeLog, AnalyzeLeeds, AnalyzePipsPro)
        for analysis_class in analysis_classes:
            analysis_instance = analysis_class(path, config)
            if analysis_instance.keyword_in_here():
                return True, analysis_instance
        return False, None


class FileAnalyzerEvent(FileSystemEventHandler):
    """Handler for file events."""

    def __init__(self, config):
        self.config = config

    def on_created(self, event):
        """Called when a file is moved into the watched directory."""
        full_file_path = osp.abspath(event.src_path)
        local_path = osp.basename(full_file_path)
        # determine if file has keyword in it
        do_analysis, analysis_instance = analysis_should_be_done(full_file_path, self.config)
        # process it if so
        if do_analysis:
            time.sleep(1)
            try:
                analysis_instance.process()
            except BaseException as e:
                logging.info(local_path + " encountered an error and was not processed." + e)
        else:
            logging.info(local_path + " was added but was not deemed a file to be analyzed.")


def start_watching(directory, config_file=None):
    """Start watching the directory and analyze any applicable files that may be moved there."""
    logging.info("Starting watcher...")
    # set up configuration
    config = load_config(config_file)
    # set up file watcher
    event_handler = FileAnalyzerEvent(config)
    observer = Observer()
    observer.schedule(event_handler, directory, recursive=True)
    logging.info("Pylinac now watching at " + osp.abspath(directory))
    observer.start()
    while True:
        time.sleep(1)
    observer.join()


def load_config(config_file=None):
    """Load an external configuration YAML file, or load the default one."""
    if config_file is None:
        yaml_config_file = osp.join(osp.dirname(__file__), 'watcher_config.yaml')
    else:
        yaml_config_file = config_file
    config = yaml.load(open(yaml_config_file).read())
    return config
