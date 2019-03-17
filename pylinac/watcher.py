"""The watcher file is a script meant to be run as an ongoing process to watch a given directory and analyzing files
that may be moved there for certain keywords. Automatic processing will be started if the file contains the keywords."""
import collections
import concurrent.futures
import datetime
import os
from functools import lru_cache
import gzip
import importlib
import logging
import os.path as osp
import pickle
import shutil
import time
import zipfile

import yagmail
import yaml

from pylinac.core.decorators import value_accept
from pylinac.core.io import retrieve_demo_file, is_dicom_image
from pylinac.core.image import prepare_for_classification, DicomImage
# from pylinac.core import schedule

from pylinac import DRMLC, DRGS, Starshot, PicketFence, WinstonLutz, LeedsTOR, StandardImagingQC3, load_log, LasVegas
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
    """
    obj = object
    config_name = ''
    has_classification = False

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

    @classmethod
    def run(cls, files, config, skip_list):
        files = drop_skips(files, skip_list)
        for file in files:
            cond1 = contains_keywords(file, config, cls.config_name)
            if config[cls.config_name]['use-classifier']:
                cond2 = matches_classifier(file, cls)
            else:
                cond2 = False
            if cond1 or cond2:
                obj = cls(file, config)
                obj.process()
                skip_list.append(osp.basename(file))

    def process(self):
        """Process the file; includes analysis, saving results to file, and sending emails."""
        logger.info(self.full_path + " will be analyzed...")
        self.instance = self.obj(self.full_path, **self.constructor_kwargs)
        self.analyze()
        if self.config['email']['enable-all']:
            self.send_email()
        elif self.config['email']['enable-failure'] and self.should_send_failure_email():
            self.send_email()
        self.publish_pdf()
        # self.save_zip()
        logger.info("Finished analysis on " + self.local_path)

    def save_zip(self):
        # save results and original file to a compressed ZIP archive
        with zipfile.ZipFile(self.zip_filename, 'w', compression=zipfile.ZIP_DEFLATED) as zfile:
            zfile.write(self.full_path, arcname=osp.basename(self.full_path))
        # remove the original files
        os.remove(self.full_path)

    @property
    def constructor_kwargs(self):
        """Any keyword arguments meant to be given to the constructor call."""
        return {}

    @property
    def zip_filename(self):
        """The name of the file for the ZIP archive."""
        return self.base_name + self.config['general']['file-suffix'] + '.zip'

    @property
    def pdf_filename(self):
        """The name of the file for the PDF results."""
        return self.base_name + '.pdf'

    @property
    def keywords(self):
        """The keywords that signal a file is of a certain analysis type."""
        return self.config[self.config_name]['keywords']

    def keyword_in_here(self):
        """Determine whether a keyword exists in the filename."""
        return any(keyword in self.local_path.lower() for keyword in self.keywords)

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
            attachments = [self.pdf_filename]
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

    def publish_pdf(self):
        self.instance.publish_pdf(self.pdf_filename, unit=self.config['general']['unit'])

    def should_send_failure_email(self):
        """Check whether analysis results were poor and an email should be triggered."""
        return not self.instance.passed

    def analyze(self):
        """Analyze the file."""
        self.instance.analyze(**self.analysis_settings)


class AnalyzeLeeds(AnalyzeMixin):
    """Analysis runner for Leeds TOR phantom."""
    obj = LeedsTOR
    config_name = 'leeds'


class AnalyzeQC3(AnalyzeMixin):
    """Analysis runner for Standard Imaging QC-3."""
    obj = StandardImagingQC3
    config_name = 'qc3'


class AnalyzeLasVegas(AnalyzeMixin):
    """Analysis runner for Las Vegas phantom."""
    obj = LasVegas
    config_name = 'las-vegas'


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


class AnalyzeLog(AnalyzeMixin):
    """Analysis runner for dynalogs or trajectory logs."""
    obj = load_log
    config_name = 'logs'

    @property
    def log_time(self):
        rev_path = self.local_path[::-1]
        u_idx = rev_path.find('_')
        log_time = osp.splitext(rev_path[:u_idx][::-1])[0]
        return log_time

    def send_email(self):
        """Send an email with the analysis results."""
        super().send_email(name=self.log_time, attachments='')

    def analyze(self):
        """Log analysis is done via calculating gamma."""
        self.instance.fluence.gamma.calc_map(**self.analysis_settings)

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
            self.publish_pdf()
            self.save_zip()
            if self.config['email']['enable-all']:
                self.send_email()
            elif self.config['email']['enable-failure'] and self.should_send_failure_email():
                self.send_email()
        logger.info("Finished analysis on " + self.local_path)
        return True

    @classmethod
    def run(cls, files, config, skip_list):
        files = drop_skips(files, skip_list)
        for file in files:
            if contains_keywords(file, config, cls.config_name):
                obj = cls(file, config)
                obj.process()
                skip_list.append(osp.basename(file))


class AnalyzeCatPhan(AnalyzeMixin):
    """Analysis runner for CBCTs."""
    config_name = 'catphan'

    def __init__(self, path, config, zip_format=True):
        self.zip_format = zip_format
        super().__init__(path=path, config=config)
        p = importlib.import_module('pylinac')
        if zip_format:
            self.obj = getattr(p, config[self.config_name]['model']).from_zip
        else:
            self.obj = getattr(p, config[self.config_name]['model'])

    @property
    def pdf_filename(self):
        """The name of the file for the PDF results."""
        if self.zip_format:
            return self.full_path.replace('.zip', '.pdf')
        else:
            return '{osp.dirname(self.instance.dicom_stack[0].path)}\CBCT - {self.instance.dicom_stack[0].date_created(format="%A, %I-%M-%S, %B %d, %Y")}.pdf'

    def analyze(self):
        self.instance.analyze(**self.analysis_settings)

    def should_send_failure_email(self):
        """Failure of CBCT depends on individual module performance."""
        send = False
        for key, val in self.failure_settings.items():
            if key == 'hu-passed' and not self.instance.ctp404.passed_hu:
                send = True
            if key == 'uniformity-passed' and not self.instance.ctp486.overall_passed:
                send = True
            if key == 'geometry-passed' and not self.instance.ctp404.passed_geometry:
                send = True
            if key == 'thickness-passed' and not self.instance.ctp404.passed_thickness:
                send = True
        return send

    @classmethod
    def run(cls, files, config, skip_list):
        files = drop_skips(files, skip_list)
        # analyze ZIP archives
        for file in files:
            cond1 = contains_keywords(file, config, cls.config_name)
            cond2 = file.endswith('.zip')
            if cond1 and cond2:
                obj = cls(file, config)
                obj.process()
                skip_list.append(osp.basename(file))
        # analyze directory groups
        done = False
        while not done:
            files = drop_skips(files, skip_list)
            if len(files) > 50:
                try:
                    obj = cls(osp.dirname(files[0]), config, zip_format=False)
                    obj.process()
                    if config[cls.config_name]['analysis']['zip_after']:
                        skip_list.append(osp.basename(obj.pdf_filename).replace('.pdf', '.zip'))
                    else:
                        for file in obj.instance.dicom_stack:
                            skip_list.append(osp.basename(file.path))
                except Exception as e:
                    print(e)
                    done = True
            else:
                done = True


class AnalyzeWL(AnalyzeMixin):
    """Analysis runner for Winston-Lutz images."""
    obj = WinstonLutz.from_zip
    config_name = 'winston-lutz'

    def __init__(self, path, config, zip_format=True):
        self.zip_format = zip_format
        self.config = config
        if zip_format:
            super().__init__(path, config)
            self.obj = WinstonLutz.from_zip
        else:
            self.full_path = path
            self.obj = WinstonLutz

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

    def process(self):
        """Process the file; includes analysis, saving results to file, and sending emails."""
        if self.zip_format:
            logger.info(self.full_path + " will be analyzed...")
        else:
            logger.info("Winston Lutz batch will be analyzed...")
        self.instance = self.obj(self.full_path, **self.constructor_kwargs)
        self.analyze()
        if self.config['email']['enable-all']:
            self.send_email()
        elif self.config['email']['enable-failure'] and self.should_send_failure_email():
            self.send_email()
        self.publish_pdf()
        if not self.zip_format:
            self.save_zip()
        logger.info("Finished Winston-Lutz analysis")

    @property
    def pdf_filename(self):
        """The name of the file for the PDF results."""
        if self.zip_format:
            return self.base_name + '.pdf'
        else:
            dirname = osp.dirname(self.full_path[0])
            dcm = DicomImage(self.full_path[0])
            name = f'Winston-Lutz - {dcm.date_created()}.pdf'
            return osp.join(dirname, name)

    @classmethod
    def run(cls, files, config, skip_list):
        files = drop_skips(files, skip_list)
        # analyze ZIP archives
        for file in files:
            cond1 = contains_keywords(file, config, cls.config_name)
            cond2 = file.endswith('.zip')
            if cond1 and cond2:
                obj = cls(file, config)
                obj.process()
                skip_list.append(osp.basename(file))
        # analyze directory groups
        if config[cls.config_name]['use-classifier']:
            wlfiles = [f for f in files if matches_classifier(f, cls)]
            if len(wlfiles) > 3:
                obj = cls(wlfiles, config, zip_format=False)
                obj.process()
                for file in wlfiles:
                    skip_list.append(osp.basename(file))


class AnalyzeVMAT(AnalyzeMixin):
    """Analysis runner for VMATs."""
    config_name = 'vmat'

    def __init__(self, path, config, zip_format=True):
        self.zip_format = zip_format
        self.config = config
        if isinstance(path, str):
            super().__init__(path=path, config=config)
        else:
            self.path = path
        if zip_format:
            self.obj = VMAT.from_zip
        else:
            self.obj = VMAT

    def process(self):
        """Process the file; includes analysis, saving results to file, and sending emails."""
        if self.zip_format:
            logger.info(self.local_path + " found and will be analyzed...")
            self.instance = self.obj(self.full_path, **self.constructor_kwargs)
        else:
            logger.info("VMAT pair found and will be analyzed...")
            self.instance = self.obj(self.path, delivery_types=('open', 'dmlc'), **self.constructor_kwargs)
        self.analyze()
        if self.config['email']['enable-all']:
            self.send_email()
        elif self.config['email']['enable-failure'] and self.should_send_failure_email():
            self.send_email()
        self.publish_pdf()
        if not self.zip_format:
            self.save_zip()
        if self.zip_format:
            logger.info("Finished analysis on " + self.local_path)
        else:
            logger.info("Finished analysis on VMAT pair")

    def analyze(self):
        """Analyze the file."""
        if self.zip_format:
            self.instance.analyze(**self.analysis_settings)
        else:
            # determine test type
            classifcations = [matches_classifier(p, self) for p in self.path]
            if 'drgs' in classifcations:
                self.test_type = 'drgs'
            elif 'mlcs' in classifcations:
                self.test_type = 'drmlc'
            else:
                raise TypeError("Cannot determine the test type of the VMAT pair")
            self.instance.analyze(test=self.test_type, **self.analysis_settings)

    @property
    def pdf_filename(self):
        """The name of the file for the PDF results."""
        if self.zip_format:
            return self.base_name + '.pdf'
        else:
            dirname = osp.dirname(self.path[0])
            dcm = DicomImage(self.path[0])
            name = f'VMAT {self.test_type.upper()} - {dcm.date_created()}.pdf'
            return osp.join(dirname, name)

    @classmethod
    def run(cls, files, config, skip_list):
        files = drop_skips(files, skip_list)
        # analyze ZIP archives
        for file in files:
            cond1 = contains_keywords(file, config, cls.config_name)
            cond2 = file.endswith('.zip')
            if cond1 and cond2:
                obj = cls(file, config)
                obj.process()
                skip_list.append(osp.basename(file))
        # analyze directory groups
        if config[cls.config_name]['use-classifier']:
            imgs = [DicomImage(f) for f in files if is_dicom_image(f)]
            uids = [i.metadata.SeriesInstanceUID for i in imgs]
            c = collections.Counter(uids)
            for uid, num in c.items():
                if num == 2:
                    img_uids = [i.path for i in imgs if i.metadata.SeriesInstanceUID == uid]
                    obj = cls(img_uids, config, zip_format=False)
                    obj.process()
                    for img in img_uids:
                        skip_list.append(osp.basename(img))


def drop_skips(files, skip_list):
    return [f for f in files if (osp.basename(f) not in skip_list)]


def matches_classifier(file, obj):
    if obj in (AnalyzePF, AnalyzeQC3, AnalyzeLeeds, AnalyzeStar, AnalyzeWL):
        img_type = 'single'
        match = {
            AnalyzePF: 1,
            AnalyzeQC3: 2,
            AnalyzeLeeds: 3,
            AnalyzeStar: 4,
            AnalyzeWL: 5,
        }
    elif isinstance(obj, AnalyzeVMAT):
        img_type = 'vmat'
        match = {
            1: 'open',
            2: 'drgs',
            3: 'mlcs',
        }
    clf = get_image_classifier(img_type)
    try:
        img = prepare_for_classification(file)
        classification = clf.predict(img.reshape(1, -1))
    except:
        return False
    else:
        if img_type == 'single':
            return classification[0] == match[obj]
        if img_type == 'vmat':
            return match[classification[0]]


def contains_keywords(file, config, obj):
    return any(keyword.lower() in osp.basename(file.lower()) for keyword in config[obj]['keywords'])


def _copy_new_files(directory, config):
    """Move any new files from the source folder(s) into the destination folder"""
    def move(sfile, source_dir, directory):
        """The function to move files from the source to the destination"""
        sbase = osp.splitext(osp.split(sfile)[1])[0]
        already_here = any(sbase in dfile for dfile in dest_files)
        if not already_here:
            shutil.copy(osp.join(source_dir, sfile), directory)
            logger.info(f"Copied {sfile} into pylinac directory")

    # return if no sources are configured or are not real directories
    source_undefined = config['general']['sources'] is None
    if source_undefined:
        return
    sources_are_dirs = all(osp.isdir(source) for source in config['general']['sources'])
    if not sources_are_dirs:
        return

    # move new files into destination directory
    dest_files = os.listdir(directory)
    with concurrent.futures.ThreadPoolExecutor(4) as exec:
        for source_dir in config['general']['sources']:
            logger.info(f"Querying new files from {source_dir}")
            source_files = os.listdir(source_dir)
            time.sleep(0.5)
            if config['general']['rolling-window-days'] > 0:
                window_cutoff = datetime.datetime.timestamp(datetime.datetime.now() - datetime.timedelta(days=config['general']['rolling-window-days']))
                source_files = [f for f in source_files if osp.getmtime(osp.join(source_dir, f)) > window_cutoff]
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


def analyze_new_files(directory, config, force):
    """Analyze any new files that have been moved into the directory pylinac is watching.

    Parameters
    ----------
    directory : str
        The directory pylinac is watching
    config : str
        The YAML configuration file; see :func:`~pylinac.watcher.load_config`.
    """
    if force:
        skip_list = []
    else:
        skip_list = get_skip_list(directory)
    for pdir, _, files in os.walk(directory):
        # drop files with avoidance keywords
        files = [osp.join(pdir, f) for f in files if not any((k in f) for k in config['general']['avoid-keywords'])]
        # drop files outside of time window
        if config['general']['rolling-window-days'] > 0:
            window_cutoff = datetime.datetime.timestamp(datetime.datetime.now() - datetime.timedelta(days=config['general']['rolling-window-days']))
            files = [f for f in files if osp.getmtime(f) > window_cutoff]
        # analyze
        for analysis in (AnalyzeLog, AnalyzeCatPhan, AnalyzeVMAT, AnalyzeWL, AnalyzePF, AnalyzeLeeds, AnalyzeStar, AnalyzeQC3):
            analysis.run(files, config, skip_list)
    # update skip list
    set_skip_list(directory, skip_list)


# def watch(directory=None, config_file=None):
#     """Start watching the directory and analyze any applicable files that may be moved there.
#
#     Parameters
#     ----------
#     directory : str, None
#         The path to the directory that pylinac will monitor for new files and keep analysis results.
#         If None, the directory will be pulled from the config file. If no path is specified either
#         by the argument or config file an error will be raised.
#     config_file : str, None
#         The path to the YAML configuration file.
#         If None (default), will load the default config file.
#     """
#     logger.info("Starting watcher...")
#     # set up configuration
#     config = load_config(config_file, verbose=True)
#     query_freq = config['general']['query-frequency']
#     logger.info(f"Querying frequency: {query_freq:1.0f}s")
#     schedule.every(query_freq).seconds.do(process, directory, config_file, True, False)
#     while True:
#         schedule.run_pending()
#         time.sleep(1)


def process(directory=None, config_file=None, copy_new_files=False, verbose=True, force=False):
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
    copy_new_files : bool
        Whether to move new logs from the sources in the config file.
    verbose : bool
        Whether to include certain logging statements.
    force : bool
        If True, will reanalyze all files within the analysis folder regardless of the pylinac skip list.
        If False, the skip list will be respected. 
        This flag may be useful if you've changed your analysis settings and want to reanalyze everything for example.
    """
    config = load_config(config_file, verbose)
    if directory is None:
        if not osp.isdir(config['general']['directory']):
            raise NotADirectoryError("No directory was passed nor did the config file contain a valid analysis directory")
        else:
            directory = config['general']['directory']
            if verbose:
                logger.info(f"Performing analysis on directory: {directory}")
    if copy_new_files:
        _copy_new_files(directory, config)
    analyze_new_files(directory, config, force=force)


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
        logger.info(f"Using configuration file: {yaml_config_file}")
    return config


@lru_cache(maxsize=1)
@value_accept(img_type=('single', 'vmat'))
def get_image_classifier(img_type):
    """Load the CBCT HU slice classifier model. If the classifier is not locally available it will be downloaded."""
    if img_type == 'single':
        classifier = 'singleimage_classifier.pkl.gz'
    elif img_type == 'vmat':
        classifier = 'vmat_classifier.pkl.gz'
    classifier_file = osp.join(osp.dirname(__file__), 'demo_files', classifier)
    # get the classifier if it's not downloaded
    if not osp.isfile(classifier_file):
        logger.info("Downloading classifier from the internet...")
        classifier_file = retrieve_demo_file(classifier)
        logger.info("Done downloading")

    with gzip.open(classifier_file, mode='rb') as m:
        clf = pickle.load(m)
    return clf
