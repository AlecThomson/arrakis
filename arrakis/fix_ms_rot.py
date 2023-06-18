#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fix the FEED rotation of ASKAP MSs
"""

from pathlib import Path

from casacore import table
import numpy as np

def get_pol_axis(ms: Path) -> float:
    """
    Get the polarization axis from the ASKAP MS
    """
    with table("%s/FEED" %(ms), readonly=True, ack=False) as tf:
        ms_feed = tf.getcol('RECEPTOR_ANGLE')
        pol_axes = -(np.degrees(ms_feed)-45.0)

    assert (ms_feed[:,0] == ms_feed[0,0]).all() & (ms_feed[:,1] == ms_feed[0,1]).all(), \
        "The RECEPTOR_ANGLE changes with time, please check the MS"
    
    pol_axis = pol_axes[0,0]
    return pol_axis