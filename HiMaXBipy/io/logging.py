#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: David Kaltenbrunner
"""
import logging


def setup_logfile(log, logname):
    handler = logging.FileHandler(filename=logname, mode='w')
    handler.setLevel(logging.INFO)
    log_state = log.handlers.copy()
    log.addHandler(handler)
    return log_state


def set_loglevel(log, level):
    log.setLevel(level)
    for handler in log.handlers:
        handler.setLevel(level)


def setup_logger(name, working_dir):
    log = logging.getLogger(name)
    log.setLevel(logging.DEBUG)
    logFormatter = logging.Formatter("%(levelname)s: %(message)s")
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    consoleHandler.setLevel(logging.WARNING)
    log.addHandler(consoleHandler)
    fileHandler = logging.FileHandler(filename=f'{working_dir}'
                                      '/logfiles/last_run.log', mode='w')
    fileHandler.setFormatter(logFormatter)
    fileHandler.setLevel(logging.INFO)
    log.addHandler(fileHandler)
    return log
