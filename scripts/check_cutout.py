#!/usr/bin/env python3
from astropy.table import Table
from spectral_cube import SpectralCube
from glob import glob
import os
import matplotlib.pyplot as plt
import astropy.units as u 
from matplotlib.patches import Ellipse
from astropy.visualization import SqrtStretch, ZScaleInterval, ImageNormalize, LogStretch, MinMaxInterval

def plot_comp(row, c='b'):
    ax = plt.gca()
    ax.plot(
        row['RA'], 
        row['Dec'], 
        marker='X', 
        color=c,
        transform=ax.get_transform('world')
    )
    e = Ellipse(
        ((row['RA']*u.deg).value, (row['Dec']*u.deg).value),
        width=(row['Maj']*u.arcsec).to(u.deg).value,
        height=(row['Min']*u.arcsec).to(u.deg).value,
        angle=row['PA']+90,
        transform=ax.get_transform('world'),
        facecolor='none',
        edgecolor=c,
        linewidth=2
    )
    ax.add_artist(e)

    ax.text(
        (row['RA']*u.deg).value, (row['Dec']*u.deg).value,
        row['Peak_flux'],
        transform=ax.get_transform('world')
    )



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
    mom0 = cube.mean(axis=0)

    fig = plt.figure()
    ax = plt.subplot(projection=mom0.wcs)
    norm = ImageNormalize(mom0.value, interval=MinMaxInterval(),
                      stretch=SqrtStretch())
    im = ax.imshow(mom0.value, norm=norm)
    ax.plot(sub['RA']*u.deg, sub['Dec']*u.deg, 'rX', transform=ax.get_transform('world'))
    ax.plot(sub_S['RA']*u.deg, sub_S['Dec']*u.deg, 'bX', transform=ax.get_transform('world'))
    ra = ax.coords['ra']
    dec = ax.coords['dec']
    ra.set_major_formatter('d.dd')
    dec.set_major_formatter('d.dd')
    fig.colorbar(im, label=f'Average flux density [{mom0.unit}]')
    for row in sub:
        plot_comp(row, c='w')
    plt.show()
    # plt.savefig('test.png')



if __name__ == "__main__":
    import argparse
    descStr = """
    Check cutouts
    """
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'cutdir',
        metavar='cutdir',
        type=str,
        help='Cutout dir')

    parser.add_argument(
        'source',
        metavar='source',
        type=str,
        help='Source name.')
    
    args = parser.parse_args()
    main(
        args.cutdir,
        args.source
    )