#!/usr/bin/env python
"""JSON utilities"""

import dataclasses
import json

import numpy as np
from astropy.io import fits

from arrakis.utils.fitsutils import head2dict


class MyEncoder(json.JSONEncoder):
    """Cutom JSON encorder.

    Parses the data stored in source_dict to JSON without
    errors.

    """

    def default(self, obj):  # pylint: disable=E0202
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.complex):
            return (obj.real, obj.imag)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, fits.Header):
            return head2dict(obj)
        elif dataclasses.is_dataclass(obj):
            return dataclasses.asdict(obj)
        else:
            return super(MyEncoder, self).default(obj)
