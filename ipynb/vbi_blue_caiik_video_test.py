import numpy as np
import matplotlib.pyplot as plt
import juanfit
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.visualization import ImageNormalize, SqrtStretch, LogStretch, ZScaleInterval
from ccdproc import ImageFileCollection
from glob import glob
import os
import cmcrameri.cm as cmcm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import AutoLocator, AutoMinorLocator, FixedLocator, FixedFormatter, LogLocator, StrMethodFormatter

def plot_colorbar(im, ax, width="3%", height="100%",loc="lower left",fontsize=14):
    clb_ax = inset_axes(ax,width=width,height=height,loc=loc,
                bbox_to_anchor=(1.02, 0., 1, 1),
                 bbox_transform=ax.transAxes,
                 borderpad=0)
    clb = plt.colorbar(im,pad = 0.05,orientation='vertical',ax=ax,cax=clb_ax)
    clb_ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    clb_ax.yaxis.get_offset_text().set_fontsize(fontsize)
    clb_ax.tick_params(labelsize=fontsize)
    return clb, clb_ax


path = "../src/DKIST/pid_1_118/BJLKB/"
vbi_file_collection = ImageFileCollection(path,glob_include="*L1.fits")
vbi_file_df = vbi_file_collection.summary.to_pandas()

savepath = "../img/DKIST/vbi_caii/"

def make_plot(path,filename,savepath, savefilename):
    with fits.open(os.path.join(path,filename)) as hdul:
        hdr = hdul[1].header
        data = hdul[1].data

    xcoord = hdr["CRVAL1"] + hdr["CDELT1"]*(np.linspace(1,hdr["NAXIS1"],hdr["NAXIS1"]) - hdr["CRPIX1"])
    ycoord = hdr["CRVAL2"] + hdr["CDELT2"]*(np.linspace(1,hdr["NAXIS2"],hdr["NAXIS2"]) - hdr["CRPIX2"])

    fig,ax = plt.subplots(figsize=(8,8),constrained_layout=True)
    norm = ImageNormalize(data,interval=ZScaleInterval(),stretch=SqrtStretch())
    im = ax.pcolormesh(xcoord, ycoord, data/hdr["TEXPOSUR"],rasterized=True,cmap=cmcm.lajolla_r,
                shading="auto")

    ax.tick_params(labelsize=18,direction="in",which="both",top=True,right=True)
    ax.tick_params(which="major",length=8,width=1.5)
    ax.tick_params(which="minor",length=6,width=1.5)
    clb, clbax = plot_colorbar(im, ax, width="5%",fontsize=18)
    ax.set_xlabel(r"SOLAR-X [arcsec]",fontsize=18)
    ax.set_ylabel(r"SOLAR-Y [arcsec]",fontsize=18)
    ax.set_title(r"DKIST/VBI-Blue Ca \textsc{ii} K",fontsize=18)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_title(r"DKIST/VBI-Blue Ca \textsc{ii} K "+hdr["DATE-AVG"][10:19],fontsize=18)

    plt.savefig(fname=os.path.join(savepath, savefilename),dpi=300,format="png",
                bbox_inches="tight")
    fig.clf()
    for ax_ in (ax,clbax):
        ax_.cla()
    plt.close(fig)

for ii, row_ in vbi_file_df.iterrows():
    make_plot(path,row_["file"],savepath,"{:03d}.png".format(ii))

    

