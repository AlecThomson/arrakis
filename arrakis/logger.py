#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Logging module for arrakis"""

import logging

# Create formatter
# formatter = logging.Formatter(
#     "SPICE: %(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s"
# )
class CustomFormatter(logging.Formatter):
    grey = "\x1b[38;20m"
    blue = "\x1b[34;20m"
    green = "\x1b[32;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format_str = "%(asctime)s.%(msecs)03d %(module)s - %(funcName)s: %(message)s"

    FORMATS = {
        logging.DEBUG: f"{blue}SPICE-%(levelname)s{reset} {format_str}",
        logging.INFO: f"{green}SPICE-%(levelname)s{reset} {format_str}",
        logging.WARNING: f"{yellow}SPICE-%(levelname)s{reset} {format_str}",
        logging.ERROR: f"{red}SPICE-%(levelname)s{reset} {format_str}",
        logging.CRITICAL: f"{bold_red}SPICE-%(levelname)s{reset} {format_str}",
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt, "%Y-%m-%d %H:%M:%S")
        return formatter.format(record)


def get_arrakis_logger(name: str='arrakis') -> logging.Logger:
    """Will construct a logger object.

    Args:
        name (str, optional): Name of the logger to attempt to use. This is ignored if in a prefect flowrun. Defaults to 'arrakis'.

    Returns:
        logging.Logger: The appropriate logger
    """
    logging.captureWarnings(True)
    logger = logging.getLogger(name)
    logger.setLevel(logging.WARNING)

    # Create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    
    # Add formatter to ch
    ch.setFormatter(CustomFormatter())
    logger.addHandler(ch)
        
    return logger

logger = get_arrakis_logger()
