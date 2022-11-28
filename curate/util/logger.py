"""
    Code to handler the log and different errors that may appear.
    Copied from Flame (https://github.com/phi-grib/flame/blob/master/flame/util/logger.py) by Manuel Pastor (manuel.pastor@upf.edu)

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 25/11/2022, 12:06 PM
"""

import appdirs
import functools
import logging
import sys

from logging.handlers import RotatingFileHandler
from pathlib import Path
from typing import Callable

def supress_log(logger: logging.Logger) -> Callable:
    """
        Decorator for suprerss logs during objects workflow
        Logs we are entering a supress log routine and
        disables the logger setting the minimum message level at
        interpreter level.

        :param logger: logger object

        :return supressor: supressor function defined inside the decorator function
        :return decorator: decorator function acting recursively
    """
    def decorator(func):
        @functools.wraps(func)
        def supressor(*args, **kwargs):
            logger.warning('Entering OBJECTS workflow. Logger will be disabled'
                           ' below error level')
            logging.disable(logging.WARNING)
            func_results = func(*args, **kwargs)
            logging.disable(logging.NOTSET)
            logger.debug('Logger enabled again!')
            return func_results
        return supressor
    return decorator


def get_log_file() -> Path:
    """ 
        Returns the log file path
        The path of the log file is given by
        appdirs.user_log_dir

        :return log_filename: file name and path of the log file
    """
    
    log_filename_path = './'
    log_filename_path = Path(log_filename_path)
    
    # creeate dir if it does not exist
    if not log_filename_path.exists():
        log_filename_path.mkdir(parents=True)

    log_filename = log_filename_path / 'datacur.log'  # append file name

    # check if exists to not erase current file
    if not log_filename.exists():
        log_filename.touch()

    return log_filename


def get_logger(name: str) -> logging.Logger:
    """
        Inits a logger and adds the handlers.
        If the logger is already created doesn't adds new handlers
        since those are set at interpreter level and already exists.

        :param name: name of the code executing the logger.

        :return logger: logger obtained after calling LOG() object functions
    """

    # create logger
    logger = logging.getLogger(name)

    # set base logger level to DEBUG but fine tu the handlers
    # for custom level
    logger.setLevel(logging.DEBUG)

    # create formatter fdor file handler (more explicit)
    file_formatter = logging.Formatter(
        '%(levelname)-8s [%(asctime)s] %(thread)d - %(name)s - %(message)s'
    )

    # formater for stream handler (less info)
    stdout_formatter = logging.Formatter(
        '%(levelname)s - %(message)s'
    )

    log_file = get_log_file()  # Create the log file
    # create console and file handler
    # if not already created
    if not logger.handlers:

        # Send DEBUG to a rotating log file
        # Limit the size to 1000000Bytes ~ 1MB 
        fh = RotatingFileHandler(log_file, maxBytes=1000000, backupCount=3)
        fh.set_name('filehandler')
        fh.setLevel('DEBUG')
        fh.setFormatter(file_formatter)
        logger.addHandler(fh)

        # send INFO to the console (stdin)
        ch = logging.StreamHandler(sys.stdout)
        ch.set_name('streamhandler')
        ch.setLevel('INFO')
        ch.setFormatter(stdout_formatter)
        logger.addHandler(ch)

        return logger
        
    # if there already handlers just return the logger
    # since its already configured
    else:
        return logger