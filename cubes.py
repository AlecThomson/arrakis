import astropy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.table import Table
from astropy.io import fits
from glob import glob
import os
from tqdm import tqdm

rootdir = '/group/askap/mcc381/RACS/'
sbid = '8585/'
field = 'RACS_test4_1.05_1049-31A/'
files = glob(rootdir+sbid+field+'image.restored.i*contcube*linmos*fits')

selavyfiles = glob(rootdir+sbid+field+'selavy*/*islands*.xml')

islands = Table.read(selavyfiles[0], format='votable')

hdulist = fits.open(files[0], mode='denywrite')
hdu = hdulist[0]
data = np.squeeze(hdu.data)
hdulist.close()

island_idx = 10
cutouts = []
#for i in tqdm(range(9)):
#    cutout = np.nanmean(data[:, \
#        islands['col_y_min'][i]:islands['col_y_max'][i], \
#            islands['col_x_min'][i]:islands['col_x_max'][i]]\
#                , axis=0)
#    cutouts.append(cutout)

#pix = 100
#cutout = np.nanmean(data[:, int(islands['col_y_cen'][island_idx])-pix:int(islands['col_y_cen'][island_idx])+pix, int(islands['col_x_cen'][island_idx])-pix:int(islands['col_x_cen'][island_idx])+pix], axis=0)

#fig, ax = plt.subplots(3, 3, sharex='col', sharey='row')
#for i in range(3):
#    for j in range(3):
#        ax[i, j].imshow(cutout[i+j])
#plt.savefig('test.png')
#hdu = 