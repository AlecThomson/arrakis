#!/usr/bin/env python3
"""Logging module for arrakis."""

from __future__ import annotations

import argparse
import io
import logging

from arrakis.utils.typing import Struct


# https://stackoverflow.com/questions/61324536/python-argparse-with-argumentdefaultshelpformatter-and-rawtexthelpformatter
class UltimateHelpFormatter(
    argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter
):
    """Combines RawTextHelpFormatter and ArgumentDefaultsHelpFormatter."""


class Formats(Struct):
    """Log formats.

    Attributes:
        debug (str): Debug log format
        info (str): Info log format
        warning (str): Warning log format
        error (str): Error log format
        critical (str): Critical log format

    """

    debug: str
    """Debug log format"""
    info: str
    """Info log format"""
    warning: str
    """Warning log format"""
    error: str
    """Error log format"""
    critical: str
    """Critical log format"""


class TqdmToLogger(io.StringIO):
    """Output stream for TQDM which will output to logger module instead of the StdOut."""

    logger = None
    level = None
    buf = ""

    def __init__(self, logger, level=None):
        """TQDM logger."""
        super().__init__()
        self.logger = logger
        self.level = level or logging.INFO

    def write(self, buf):
        """Write to the buffer."""
        self.buf = buf.strip("\r\n\t ")

    def flush(self):
        """Flush the buffer."""
        self.logger.log(self.level, self.buf)


# Create formatter
# formatter = logging.Formatter(
#     "SPICE: %(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s"
# )
class CustomFormatter(logging.Formatter):
    """Custom formatter for logging."""

    grey = "\x1b[38;20m"
    blue = "\x1b[34;20m"
    green = "\x1b[32;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format_str = "%(asctime)s.%(msecs)03d %(module)s - %(funcName)s: %(message)s"

    FORMATS = Formats(
        debug=f"{blue}SPICE-%(levelname)s{reset} {format_str}",
        info=f"{green}SPICE-%(levelname)s{reset} {format_str}",
        warning=f"{yellow}SPICE-%(levelname)s{reset} {format_str}",
        error=f"{red}SPICE-%(levelname)s{reset} {format_str}",
        critical=f"{bold_red}SPICE-%(levelname)s{reset} {format_str}",
    )

    def format(self, record) -> str:
        """Format the log record.

        Args:
            record (LogRecord): The log record.

        Returns:
            str: Formatted log.
        """
        log_fmt = self.FORMATS._asdict().get(record.levelname.lower())
        formatter = logging.Formatter(log_fmt, "%Y-%m-%d %H:%M:%S")
        return formatter.format(record)


def get_arrakis_logger(
    name: str = "arrakis", attach_handler: bool = True
) -> logging.Logger:
    """Will construct a logger object.

    Args:
        name (str, optional): Name of the logger to attempt to use. This is ignored if in a prefect flowrun. Defaults to 'arrakis'.
        attach_handler (bool, optional): Attacjes a custom StreamHandler. Defaults to True.

    Returns:
        logging.Logger: The appropriate logger
    """
    logging.captureWarnings(True)
    logger = logging.getLogger(name)
    logger.setLevel(logging.WARNING)

    if attach_handler:
        # Create console handler and set level to debug
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)

        # Add formatter to ch
        ch.setFormatter(CustomFormatter())
        logger.addHandler(ch)

    return logger


logger = get_arrakis_logger()
