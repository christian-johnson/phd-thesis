import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from scipy.signal import convolve2d


filenames = ['fermi5','fermi1','fermi2','fermi3','fermi4']
for fitsfile in filenames:
    image_data = fits.getdata(fitsfile+'.fits')
    filename = get_pkg_data_filename(fitsfile+'.fits')
    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)
    fig = plt.figure(figsize=[10,7])
    ax = fig.add_subplot(111, projection=wcs)
    print(np.max(image_data))
    mappable=plt.imshow(np.sqrt(np.sqrt((image_data))),cmap='inferno',origin='lower')
    plt.savefig(fitsfile+'.pdf')
