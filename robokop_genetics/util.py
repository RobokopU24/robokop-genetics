import logging
import os
from logging.handlers import RotatingFileHandler


def get_logging_path():
    if 'ROBO_GENETICS_LOGS' in os.environ:
        return f'{os.environ["ROBO_GENETICS_LOGS"]}/'
    elif 'DATA_SERVICES_LOGS' in os.environ:
        return f'{os.environ["DATA_SERVICES_LOGS"]}/'
    elif 'ORION_LOGS' in os.environ:
        return f'{os.environ["ORION_LOGS"]}/'
    elif 'RAGS_HOME' in os.environ:
        return f'{os.environ["RAGS_HOME"]}/logs/'
    elif 'ROBOKOP_HOME' in os.environ:
        return f'{os.environ["ROBOKOP_HOME"]}/logs/'
    else:
        return None

def get_logger(name):
    """
        Logging utility controlling format and setting initial logging level
    """

    # get the logger with the specified name
    logger = logging.getLogger(name)

    # if it already has handlers, it was already instantiated - return it
    if logger.hasHandlers():
        return logger

    logger.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)-15s - %(funcName)s(): %(message)s')

    logs_dir_path = get_logging_path()
    if logs_dir_path is not None:
        # create a rotating file handler, 100mb max per file with a max number of 10 files
        file_handler = RotatingFileHandler(filename=os.path.join(logs_dir_path, name + '.log'),
                                           maxBytes=100000000, backupCount=10)

        # set the formatter
        file_handler.setFormatter(formatter)

        # add the handler to the logger
        logger.addHandler(file_handler)

    # create a stream handler as well (default to console/stdout)
    stream_handler = logging.StreamHandler()

    # set the formatter on the console stream
    stream_handler.setFormatter(formatter)

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
