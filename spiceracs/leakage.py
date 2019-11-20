import astropy
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os

rootdir = '/group/askap/mcc381/RACS/'
files = glob(rootdir+'*[0-9]/'+'RACS*/'+'image*v*beam*taylor.0*.fits')
print(files)
#field_dirs = glob(sbid_dirs+'RACS_test*/')

plt.figure()
plt.show()