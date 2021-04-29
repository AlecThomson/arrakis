from astropy.table import Table
from spectral_cube import SpectralCube
from glob import glob
import os
import matplotlib.pyplot as plt
import astropy.units as u 
from matplotlib.patches import Ellipse

def plot_comp(row, c='b'):
    ax = plt.gca()
    ax.plot(
        row['RA']*u.deg, 
        row['Dec']*u.deg, 
        marker='X', 
        color=c,
        transform=ax.get_transform('world')
    )
    e = Ellipse(
        (row['RA']*u.deg, row['Dec']*u.deg),
        width=row['Maj'],
        height=row['Min'],
        angle=row['PA'],
        transform=ax.get_transform('world')
    )
    ax.add_artist(e)



def main(cutdir, source):
    cat = Table.read(
        '/home/athomson/athomson/projects/RACS/Catalogues/RACS-25asec-Mosaiced_Gaussians_Final_GalCut_v2020_11_15.fits'
    )
    cat_S = Table.read(
        '/home/athomson/athomson/projects/RACS/Catalogues/RACS-25asec-Mosaiced_Sources_Final_GalCut_v2020_11_15.fits'
    )
    cat.add_index('Source_ID')
    cat_S.add_index('Source_ID')
    sub = Table(cat.loc[source])
    cat_S.add_index('Source_ID')
    sub_S = Table(cat_S.loc[source])

    cutdir = os.path.abspath(cutdir)
    images = glob(f"{cutdir}/{source}/*.image.*.i.*.linmos.fits")

    cube = SpectralCube.read(images[0])
    mom0 = cube.sum(axis=0)

    fig = plt.figure()
    ax = plt.subplot(projection=mom0.wcs)
    ax.imshow(mom0.value)
    ax.plot(sub['RA']*u.deg, sub['Dec']*u.deg, 'rX')
    ax.plot(sub_S['RA']*u.deg, sub_S['Dec']*u.deg, 'bX')
    ra = ax.coords['ra'] 
    dec = ax.coords['dec']
    ra.set_major_formatter('d.ddd')
    dec.set_major_formatter('d.ddd')
    plt.savefig('test.png')



if __name__ == "__main__":
    main()