import logging
import os
from logging.handlers import RotatingFileHandler


class LoggingUtil(object):
    """ Logging utility controlling format and setting initial logging level """

    @staticmethod
    def get_logging_path():
        if 'ROBO_GENETICS_HOME' in os.environ:
            return f'{os.environ["ROBO_GENETICS_HOME"]}/'
        elif 'ROBOKOP_HOME' in os.environ:
            return f'{os.environ["ROBOKOP_HOME"]}/logs/'
        elif 'RAGS_HOME' in os.environ:
            return f'{os.environ["RAGS_HOME"]}/logs/'
        else:
            return None

    @staticmethod
    def init_logging(name, level=logging.INFO, line_format='short', log_file_path=None, log_file_level=None):
        """
            Logging utility controlling format and setting initial logging level
        """
        # get a new logger
        logger = logging.getLogger(__name__)

        # is this the root
        if not logger.parent.name == 'root':
            return logger

        # define the various output formats
        format_type = {
            "short": '%(funcName)s(): %(message)s',
            "medium": '%(asctime)-15s - %(funcName)s(): %(message)s',
            "long": '%(asctime)-15s  - %(filename)s %(funcName)s() %(levelname)s: %(message)s'
        }[line_format]

        # create a stream handler (default to console)
        stream_handler = logging.StreamHandler()

        # create a formatter
        formatter = logging.Formatter(format_type)

        # set the formatter on the console stream
        stream_handler.setFormatter(formatter)

        # get the name of this logger
        logger = logging.getLogger(name)

        # set the logging level
        logger.setLevel(level)

        # if there was a file path passed in use it
        if log_file_path is not None:
            # create a rotating file handler, 100mb max per file with a max number of 10 files
            file_handler = RotatingFileHandler(filename=os.path.join(log_file_path, name + '.log'), maxBytes=100000000,
                                               backupCount=10)

            # set the formatter
            file_handler.setFormatter(formatter)

            # if a log level for the file was passed in use it
            if log_file_level is not None:
                level = log_file_level

            # set the log level
            file_handler.setLevel(level)

            # add the handler to the logger
            logger.addHandler(file_handler)

        # add the console handler to the logger
        logger.addHandler(stream_handler)

        # return to the caller
        return logger


class Text:
    """ Utilities for processing text. """

    @staticmethod
    def get_curie(text):
        return text.upper().split(':', 1)[0] if ':' in text else None

    @staticmethod
    def un_curie(text):
        return ':'.join(text.split(':', 1)[1:]) if ':' in text else text

    @staticmethod
    def short(obj, limit=80):
        text = str(obj) if obj else None
        return (text[:min(len(text), limit)] + ('...' if len(text) > limit else '')) if text else None

    @staticmethod
    def path_last(text):
        return text.split('/')[-1:][0] if '/' in text else text

    @staticmethod
    def snakify(text):
        decomma = '_'.join(text.split(','))
        dedash = '_'.join(decomma.split('-'))
        resu = '_'.join(dedash.split())
        return resu

    @staticmethod
    def upper_curie(text):
        if ':' not in text:
            return text
        p = text.split(':', 1)
        return f'{p[0].upper()}:{p[1]}'

    @staticmethod
    def get_curies_by_prefix(desired_curie: str, provided_synonyms: set):
        curies = set()
        for syn in provided_synonyms:
            if Text.get_curie(syn) == desired_curie:
                curies.add(syn)
        return curies
