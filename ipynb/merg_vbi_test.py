import numpy as np
import matplotlib.pyplot as plt
import juanfit
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.visualization import ImageNormalize, SqrtStretch, LogStretch, ZScaleInterval
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, SpectralCoord
from astropy.wcs.utils import celestial_frame_to_wcs
import astropy.units as u
from sunpy.coordinates.frames import Helioprojective
from ccdproc import ImageFileCollection
from glob import glob
import os
import cmcrameri.cm as cmcm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import AutoLocator, AutoMinorLocator, FixedLocator, FixedFormatter, LogLocator, StrMethodFormatter
import dkist
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
from reproject import reproject_exact, reproject_exact, reproject_interp
import copy
import h5py

def make_vbi_mosaic(tiled_ds,index):

    target_wcs = copy.deepcopy(WCS(tiled_ds[0,0].headers[0]).celestial)
    target_wcs.wcs.crpix = 6144.5, 6144.5
    target_wcs.wcs.crval = -619.0/3600., -413.0/3600.

    wcs_list = []
    for ii in range(3):
        for jj in range(3):
            if (ii == 3) and (jj ==3):
                pass
            else:
                data = tiled_ds[ii,jj][index,:,:].data.compute()
                # wcs = tiled_ds[ii,jj][index,:,:].wcs
                wcs = WCS(tiled_ds[ii,jj].headers[index]).celestial
                wcs_list.append((data,wcs))
    
    target_shape = (3*4096, 3*4096)

    array, footprint = reproject_and_coadd(wcs_list,
                                           target_wcs, shape_out=target_shape,
                                           reproject_function=reproject_interp,
                                           match_background=False)
    
    return target_wcs, target_shape, array, footprint


# asdf_file = "../src/DKIST/pid_1_118/BJLKB/VBI_L1_20220602T172155_BJLKB.asdf"
# asdf_file = "../src/DKIST/pid_1_118/BKJYA/VBI_L1_20220602T172222_BKJYA.asdf"
asdf_file = "../src/DKIST/pid_1_118/BLKGA/VBI_L1_20220602T172250_BLKGA.asdf"

ds = dkist.Dataset.from_asdf(asdf_file)

target_wcs, target_shape, vbi_merg, footprint = make_vbi_mosaic(ds,5)
# plt.imshow(vbi_merg,vmin=0,vmax=2.5e4)
# plt.show()

# with h5py.File("../sav/VBI_CaII_merg_test.h5","w") as hf:
# with h5py.File("../sav/VBI_Gband_test.h5","w") as hf:
with h5py.File("../sav/Hbeta_test.h5","w") as hf:
    df_vbi = hf.create_dataset("vbi_merg", data=vbi_merg)
    df_footprint = hf.create_dataset("footprint", data=footprint)





