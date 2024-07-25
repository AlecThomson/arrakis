"""JSON utilities."""

from __future__ import annotations

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
        """Custom JSON encoder.

        Args:
            obj (Any): Object to encode.

        Returns:
            Any: Encoded object.

        """
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, complex):
            return (obj.real, obj.imag)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, fits.Header):
            return head2dict(obj)
        if dataclasses.is_dataclass(obj):
            return dataclasses.asdict(obj)

        return super().default(obj)
