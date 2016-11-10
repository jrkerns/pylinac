"""The watcher file is a script meant to be run as an ongoing process to watch a given directory and analyzing files
that may be moved there for certain keywords. Automatic processing will be started if the file contains the keywords."""
import concurrent.futures
import datetime
import os
from functools import lru_cache
import gc
import gzip
import logging
import os.path as osp
import pickle
import shutil
import time
import zipfile

import yagmail
import yaml

from pylinac.core.io import retrieve_demo_file, retrieve_filenames
from pylinac.core.image import prepare_for_classification
from pylinac.core import schedule

from pylinac import CBCT, VMAT, Starshot, PicketFence, WinstonLutz, LeedsTOR, StandardImagingQC3, load_log
from pylinac.log_analyzer import IMAGING

logger = logging.getLogger("pylinac")

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
    save_text_method : str
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
        self.base_name = osp.splitext(self.full_path)[0]
        self.config = config

    def process(self):
        """Process the file; includes analysis, saving results to file, and sending emails."""
        logger.info(self.local_path + " will be analyzed...")
        self.instance = self.obj(self.full_path, **self.constructor_kwargs)
        self.analyze()
        if self.config['email']['enable-all']:
            self.send_email()
        elif self.config['email']['enable-failure'] and self.should_send_failure_email():
            self.send_email()
        self.save_zip()
        logger.info("Finished analysis on " + self.local_path)

    def save_zip(self):
        # save the image and/or text to file
        self.save_image()
        self.save_text()
        # save results and original file to a compressed ZIP archive
        with zipfile.ZipFile(self.zip_filename, 'w', compression=zipfile.ZIP_DEFLATED) as zfile:
            zfile.write(self.img_filename, arcname=osp.basename(self.img_filename))
            try:
                zfile.write(self.txt_filename, arcname=osp.basename(self.txt_filename))
            except:
                pass
            zfile.write(self.full_path, arcname=osp.basename(self.full_path))
        # remove the original files
        for file in (self.img_filename, self.txt_filename, self.full_path):
            try:
                os.remove(file)
            except:
                pass

    @property
    def constructor_kwargs(self):
        """Any keyword arguments meant to be given to the constructor call."""
        return {}

    @property
    def zip_filename(self):
        """The name of the file for the ZIP archive."""
        return self.base_name + self.config['general']['file-suffix'] + '.zip'

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

    def send_email(self, name=None, attachments=None):
        """Send an email with the analysis results."""
        if name is None:
            name = self.local_path
        if attachments is None:
            attachments = [self.img_filename, self.txt_filename]
        elif attachments == '':
            attachments = []
        # compose message
        current_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        if self.config['email']['enable-all']:
            statement = 'The pylinac watcher analyzed the file named or containing "{}" at {}. '
            statement += 'The analysis results are in the folder "{}".'
        elif self.config['email']['enable-failure']:
            statement = 'The pylinac watcher analyzed the file named or containing "{}" at {} and '
            statement += 'found something that failed your configuration settings.'
            statement += 'The analysis results are in the folder "{}".'
        statement = statement.format(name, current_time, osp.dirname(self.full_path))
        # send the email
        contents = [statement] + attachments
        recipients = [recipient for recipient in self.config['email']['recipients']]
        yagserver = yagmail.SMTP(self.config['email']['sender'], self.config['email']['sender-password'])
        yagserver.send(to=recipients,
                       subject=self.config['email']['subject'],
                       contents=contents)
        logger.info("An email was sent to the recipients with the results")

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
            if key == 'gantry-iso-size' and (self.instance.gantry_iso_size > val):
                send = True
            if key == 'mean-cax-bb-distance' and (self.instance.cax2bb_distance() > val):
                send = True
            if key == 'max-cax-bb-distance' and (self.instance.cax2bb_distance('max') > val):
                send = True
        return send


class AnalyzeLeeds(AnalyzeMixin):
    """Analysis runner for Leeds TOR phantom."""
    obj = LeedsTOR
    config_name = 'leeds'

    def save_text(self):
        """LeedsTOR analysis does not have a text file result"""
        pass


class AnalyzePipsPro(AnalyzeMixin):
    """Analysis runner for PipsPro QC-3."""
    obj = StandardImagingQC3
    config_name = 'pipspro'

    def save_text(self):
        """Pipspro analysis does not have a text file result"""
        pass


class AnalyzeStar(AnalyzeMixin):
    """Analysis runner for starshots."""
    obj = Starshot
    config_name = 'starshot'

    @property
    def constructor_kwargs(self):
        """Give the SID to the algorithm"""
        return {'sid': self.config[self.config_name]['analysis']['sid']}

    @property
    def analysis_settings(self):
        """Starshot analysis settings"""
        settings = super().analysis_settings
        settings.pop('sid')
        return settings


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
    obj = load_log
    save_image_method = 'save_summary'
    config_name = 'logs'

    @property
    def log_time(self):
        rev_path = self.local_path[::-1]
        u_idx = rev_path.find('_')
        log_time = osp.splitext(rev_path[:u_idx][::-1])[0]
        return log_time

    @property
    def log_txt_filename(self):
        return self.instance.txt_filename

    def send_email(self):
        """Send an email with the analysis results."""
        super().send_email(name=self.log_time, attachments='')

    def analyze(self):
        """Log analysis is done via calculating gamma."""
        self.instance.fluence.gamma.calc_map(**self.analysis_settings)

    def save_text(self):
        """Special text save method."""
        with open(self.txt_filename, 'w') as txtfile:
            txtfile.write(self.instance.report_basic_parameters(printout=False))

    def should_send_failure_email(self):
        """Failure is based on several varying criteria."""
        send = False
        for key, val in self.failure_settings.items():
            if key == 'gamma':
                gamma_below_threshold = self.instance.fluence.gamma.pass_prcnt < val
                if gamma_below_threshold:
                    send = True
            elif key == 'avg-rms':
                rms_above_threshold = self.instance.axis_data.mlc.get_RMS_avg() > val
                if rms_above_threshold:
                    send = True
            elif key == 'max-rms':
                rms_above_threshold = self.instance.axis_data.mlc.get_RMS_max() > val
                if rms_above_threshold:
                    send = True
        return send

    def process(self):
        """Process the file; includes analysis, saving results to file, and sending emails."""
        logger.info(self.local_path + " will be analyzed...")
        self.instance = load_log(self.full_path, **self.constructor_kwargs)
        if self.instance.treatment_type == IMAGING:
            logger.info(self.local_path + " is an imaging log...")
        else:
            self.analyze()
            self.save_zip()
            if self.config['email']['enable-all']:
                self.send_email()
            elif self.config['email']['enable-failure'] and self.should_send_failure_email():
                self.send_email()
        logger.info("Finished analysis on " + self.local_path)
        return True

    def save_zip(self):
        """Save the log and resulting PNG & text file to a ZIP archive."""
        # save the image and/or text to file
        self.save_image()
        self.save_text()
        # save results and original file to a compressed ZIP archive
        with zipfile.ZipFile(self.zip_filename, 'w', compression=zipfile.ZIP_DEFLATED) as zfile:
            for file in (self.img_filename, self.txt_filename, self.full_path):
                zfile.write(file, arcname=osp.basename(file))
            if self.log_txt_filename is not None:
                zfile.write(self.log_txt_filename, osp.basename(self.log_txt_filename))
        # remove the original files
        for file in (self.img_filename, self.txt_filename, self.full_path):
            os.remove(file)
        if self.log_txt_filename is not None:
            os.remove(self.log_txt_filename)


def analysis_should_be_done(path, config, skip_list):
    """Return boolean of whether the file should be analysed, based on if the filename has a keyword."""
    do_analysis, clfy = False, None
    # return if not an analysis-worthy file
    if osp.basename(path) in skip_list:
        return False, None
    has_avoid_keyword = any(item in path for item in config['general']['avoid-keywords'])
    has_suffix = config['general']['file-suffix'] in path
    if has_avoid_keyword or has_suffix:
        return False, None
    if config['general']['use-classifier']:
        do_analysis, clfy = auto_classify(path, config)
        if clfy:
            logger.info("{} classified using SVM".format(osp.basename(path)))
    if not config['general']['use-classifier'] or (clfy is None):
        do_analysis, clfy = filename_classify(path, config)
        if clfy:
            logger.info("{} classified using filename keywords".format(osp.basename(path)))
    return do_analysis, clfy


def move_logs(directory, config):
    """Move any new files from the TrueBeam log folders into the destination folder"""
    def move(sfile, source_dir, directory):
        """The function to move files from the source to the destination"""
        sbase = osp.splitext(osp.split(sfile)[1])[0]
        already_here = any(sbase in dfile for dfile in dest_files)
        if not already_here and (sfile.endswith('.bin') or sfile.endswith('.txt')):
            shutil.copy(osp.join(source_dir, sfile), directory)
            logger.info("Copied {} into pylinac directory".format(sfile))

    # return if no sources are configured or are not real directories
    log_source_undefined = config['logs']['sources'] is None
    if log_source_undefined:
        return
    log_sources_are_dirs = all(osp.isdir(source) for source in config['logs']['sources'])
    if not log_sources_are_dirs:
        return

    # move new files into destination directory
    dest_files = os.listdir(directory)
    with concurrent.futures.ThreadPoolExecutor(4) as exec:
        for source_dir in config['logs']['sources']:
            logger.info("Querying new logs from {}".format(source_dir))
            source_files = os.listdir(source_dir)
            time.sleep(0.5)
            for file in source_files:
                exec.submit(move, file, source_dir, directory)


def get_skip_list(directory):
    skip_list_file = osp.join(directory, 'pylinac.pkl')
    if osp.isfile(skip_list_file):
        with open(skip_list_file, 'rb') as sl:
            skip_list = pickle.load(sl)
    else:
        skip_list = []
    return skip_list


def set_skip_list(directory, skip_list):
    skip_list_file = osp.join(directory, 'pylinac.pkl')
    with open(skip_list_file, 'wb') as sl:
        pickle.dump(skip_list, sl)


def analyze_new_files(directory, config):
    """Analyze any new files that have been moved into the directory pylinac is watching.

    Parameters
    ----------
    directory : str
        The directory pylinac is watching
    config : str
        The YAML configuration file; see :func:`~pylinac.watcher.load_config`.
    """
    # get files that should be analyzed
    all_files = retrieve_filenames(directory)
    skip_list = get_skip_list(directory)

    def analyze(file, config, skip_list):
        """The actual analysis and processing function"""
        add2skip = False
        # determine if file has keyword in it
        do_analysis, analysis_instance = analysis_should_be_done(file, config, skip_list)
        # process it if so
        if do_analysis:
            time.sleep(0.005)
            try:
                add2skip = analysis_instance.process()
            except BaseException as e:
                logger.info(file + " encountered an error and was not processed: {}".format(e))
            gc.collect()
        return add2skip

    # start a thread pool and execute the analysis; thread pool set to 1 due to strange errors w/ >1
    for file in all_files:
        add2skip = analyze(file, config, skip_list)
        if add2skip:
            skip_list.append(osp.basename(file))

    # update skip list
    set_skip_list(directory, skip_list)


def watch(directory=None, config_file=None):
    """Start watching the directory and analyze any applicable files that may be moved there.

    Parameters
    ----------
    directory : str, None
        The path to the directory that pylinac will monitor for new files and keep analysis results.
        If None, the directory will be pulled from the config file. If no path is specified either
        by the argument or config file an error will be raised.
    config_file : str, None
        The path to the YAML configuration file.
        If None (default), will load the default config file.
    """
    logger.info("Starting watcher...")
    # set up configuration
    config = load_config(config_file, verbose=True)
    query_freq = config['general']['query-frequency']
    logger.info("Querying frequency: {:1.0f}s".format(query_freq))
    schedule.every(query_freq).seconds.do(process, directory, config_file, True, False)
    while True:
        schedule.run_pending()
        time.sleep(1)


def process(directory=None, config_file=None, move_new_logs=False, verbose=True):
    """Process the contents of the directory once through. This is contrasted to
    :func:`~pylinac.watcher.watch` which continually watches the directory.

    Parameters
    ----------
    directory : str
        The path to the directory that pylinac will process.
        If None, the directory will be pulled from the config file. If no path is specified either
        by the argument or config file an error will be raised.
    config_file : str, None
        The path to the YAML configuration file.
        If None (default), will load the default config file.
    move_new_logs : bool
        Whether to move new logs from the sources in the config file.
    verbose : bool
        Whether to include certain logging statements.
    """
    config = load_config(config_file, verbose)
    if directory is None:
        if not osp.isdir(config['general']['directory']):
            raise NotADirectoryError("No directory was passed nor did the config file contain a valid analysis directory")
        else:
            directory = config['general']['directory']
            if verbose:
                logger.info("Performing analysis on directory: {}".format(directory))
    if move_new_logs:
        move_logs(directory, config)
    analyze_new_files(directory, config)


def load_config(config_file=None, verbose=False):
    """Load a configuration YAML file.

    Parameters
    ----------
    config_file : str, None
        The path to the YAML configuration file.
        If None (default), will load the default config file.
    verbose : bool
        Whether to print more logging statements.
    """
    if config_file is None:
        yaml_config_file = osp.join(osp.dirname(__file__), 'watcher_config.yml')
    else:
        yaml_config_file = config_file
    config = yaml.load(open(yaml_config_file).read())
    if verbose:
        logger.info("Using configuration file: {}".format(yaml_config_file))
    return config


def auto_classify(path, config):
    """Classify an image using an SVM classifier.

    Parameters
    ----------
    path : str
        The path to the file to be classified.
    config : str
        The configuration settings.
    """
    clf = get_image_classifier()
    try:
        img = prepare_for_classification(path)
        classification = clf.predict(img.reshape(1, -1))
    except:
        return False, None
    else:
        if classification == 1:
            return True, AnalyzePF(path, config)
        elif classification == 2:
            return True, AnalyzePipsPro(path, config)
        elif classification == 3:
            return True, AnalyzeLeeds(path, config)
        elif classification == 4:
            return True, AnalyzeStar(path, config)
        else:
            return False, None


@lru_cache(maxsize=1)
def get_image_classifier():
    """Load the CBCT HU slice classifier model. If the classifier is not locally available it will be downloaded."""
    classifier_file = osp.join(osp.dirname(__file__), 'demo_files', 'singleimage_classifier.pkl.gz')
    # get the classifier if it's not downloaded
    if not osp.isfile(classifier_file):
        logger.info("Downloading classifier from the internet...")
        classifier_file = retrieve_demo_file('singleimage_classifier.pkl.gz')
        logger.info("Done downloading")

    with gzip.open(classifier_file, mode='rb') as m:
        clf = pickle.load(m)
    return clf


def filename_classify(path, config):
    """Classify an image using the filename convention."""
    analysis_classes = (AnalyzeStar, AnalyzeCBCT, AnalyzeVMAT, AnalyzePF, AnalyzeWL, AnalyzeLog, AnalyzeLeeds, AnalyzePipsPro)
    for analysis_class in analysis_classes:
        analysis_instance = analysis_class(path, config)
        if analysis_instance.keyword_in_here():
            return True, analysis_instance
    return False, None
