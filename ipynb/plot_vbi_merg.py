import numpy as np
import matplotlib
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
import sunpy.visualization.colormaps as cm
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

# cm_sdoaia171 = matplotlib.colormaps['sdoaia171']
# cm_eui174 = matplotlib.colormaps['solar orbiterhri_euv174']

CaII_asdf_file = "../src/DKIST/pid_1_118/BJLKB/VBI_L1_20220602T172155_BJLKB.asdf"
Gband_asdf_file = "../src/DKIST/pid_1_118/BKJYA/VBI_L1_20220602T172222_BKJYA.asdf"
Hbeta_asdf_file = "../src/DKIST/pid_1_118/BLKGA/VBI_L1_20220602T172250_BLKGA.asdf"

CaII_merg_file = "../sav/VBI_CaII_merg_test.h5"
Gband_merg_file = "../sav/VBI_Gband_test.h5"
Hbeta_merg_file = "../sav/Hbeta_test.h5"

with h5py.File(CaII_merg_file,"r") as hf:
    CaII_image = hf["vbi_merg"][:]
with h5py.File(Gband_merg_file,"r") as hf:
    Gband_image = hf["vbi_merg"][:]
with h5py.File(Hbeta_merg_file,"r") as hf:
    Hbeta_image = hf["vbi_merg"][:]

ds_CaII = dkist.Dataset.from_asdf(CaII_asdf_file)
ds_Gband = dkist.Dataset.from_asdf(Gband_asdf_file)
ds_Hbeta = dkist.Dataset.from_asdf(Hbeta_asdf_file)

wcs_CaII = copy.deepcopy(WCS(ds_CaII[0,0].headers[0]).celestial)
wcs_CaII.wcs.crpix = 6144.5, 6144.5
wcs_CaII.wcs.crval = -619.0/3600., -413.0/3600.

wcs_Gband = copy.deepcopy(WCS(ds_Gband[0,0].headers[0]).celestial)
wcs_Gband.wcs.crpix = 6144.5, 6144.5
wcs_Gband.wcs.crval = -619.0/3600., -413.0/3600.

wcs_Hbeta = copy.deepcopy(WCS(ds_Hbeta[0,0].headers[0]).celestial)
wcs_Hbeta.wcs.crpix = 6144.5, 6144.5
wcs_Hbeta.wcs.crval = -619.0/3600., -413.0/3600.

fig = plt.figure(figsize=(12,5),constrained_layout=True)

# norm_Gband = ImageNormalize(Gband_image,vmin=2e4,vmax=6e4,stretch=SqrtStretch())
ax1 = fig.add_subplot(1,3,1,projection=wcs_Gband)
ax1.imshow(Gband_image, origin="lower", cmap="Greys_r",vmin=1.5e4,vmax=5e4)
ax1.set_title(r"\textbf{G-band "+ds_Gband[0,0].headers[5]["DATE-AVG"][10:19] + r"}",fontsize=18)

ax2 = fig.add_subplot(1,3,2,projection=wcs_Hbeta)
ax2.imshow(Hbeta_image, origin="lower", cmap="Greys_r",vmin=1e4,vmax=6e4)
ax2.set_title(r"\textbf{H}$\boldsymbol{\beta}$ \textbf{"+ds_Hbeta[0,0].headers[5]["DATE-AVG"][10:19] + r"}",fontsize=18)

ax3 = fig.add_subplot(1,3,3,projection=wcs_CaII)
ax3.imshow(CaII_image, origin="lower", cmap="Greys_r",vmin=0,vmax=2.2e4)
ax3.set_title(r"\textbf{Ca \textsc{ii} K"+ds_CaII[0,0].headers[5]["DATE-AVG"][10:19] + r"}",fontsize=18)

for ax_ in (ax1, ax2, ax3):
    ax_.tick_params(labelsize=16,direction="in")
    ax_.set_ylabel(" ")
    ax_.set_xlabel(r"\textbf{SOLAR-X [arcsec]}",fontsize=16)
    ax_.grid("on")
ax1.set_ylabel(r"\textbf{SOLAR-Y [arcsec]}",fontsize=16)

plt.savefig(fname="../img/DKIST/vbi_mosaic_test.png",dpi=1000,format="png",bbox_inches="tight")