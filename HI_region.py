from __future__ import division
import numpy as np
import cPickle
import matplotlib
matplotlib.rcParams['backend'] = "Qt4Agg"
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import wcs
import copy

def make_wcs(wcs_fn):
    #Set wcs transformation
    w = wcs.WCS(wcs_fn, naxis=2)
    
    return w

def xy_to_radec(x, y, w):
    
    #Transformation
    xy = [[x, y]]
    radec = w.wcs_pix2world(xy, 1)
    
    ra = radec[0,0]
    dec = radec[0,1]
    
    return ra, dec

def radec_to_xy(ra, dec, w):
    
    #Transformation
    radec = [[ra, dec]]
    xy = w.wcs_world2pix(radec, 1)
    
    x = xy[0,0]
    y = xy[0,1]
    
    return x, y

def xcutout_data(big_galfa_data, big_galfa_hdr, xstart = 0, xstop = None):
    
    cutout_galfa_hdr = copy.copy(big_galfa_hdr)
    
    if xstop is None:
        cutout_galfa_data = big_galfa_data[:, xstart:]
    else:
        cutout_galfa_data = big_galfa_data[:, xstart:xstop]
    
    ny, nx = cutout_galfa_data.shape
    
    # Need to change the center pixels
    old_crpix1 = big_galfa_hdr["CRPIX1"]
    old_crpix2 = big_galfa_hdr["CRPIX2"]

    cutout_crpix1 = old_crpix1 - xstart

    # Define cutout header
    cutout_galfa_hdr["NAXIS"] = 2
    cutout_galfa_hdr["NAXIS1"] = nx
    cutout_galfa_hdr["NAXIS2"] = ny
    cutout_galfa_hdr["CRPIX1"] = cutout_crpix1
    
    return cutout_galfa_hdr, cutout_galfa_data
    

def ycutout_data(big_galfa_data, big_galfa_hdr, ystart = 0, ystop = None):
    
    cutout_galfa_hdr = copy.copy(big_galfa_hdr)
    
    if ystop is None:
        cutout_galfa_data = big_galfa_data[ystart:, :]
    else:
        cutout_galfa_data = big_galfa_data[ystart:ystop, :]
    
    ny, nx = cutout_galfa_data.shape
    
    # Need to change the center pixels
    old_crpix2 = big_galfa_hdr["CRPIX2"]

    cutout_crpix2 = old_crpix2 - ystart

    # Define cutout header
    cutout_galfa_hdr["NAXIS"] = 2
    cutout_galfa_hdr["NAXIS1"] = nx
    cutout_galfa_hdr["NAXIS2"] = ny
    cutout_galfa_hdr["CRPIX2"] = cutout_crpix2
    
    return cutout_galfa_hdr, cutout_galfa_data
    
nhidata_fn = "/Volumes/DataDavy/GALFA/DR2/NHIMaps/GALFA-HI_VLSR-036+0037kms_NHImap_noTcut.fits"
nhi_hdr = fits.getheader(nhidata_fn)
nhi_data = fits.getdata(nhidata_fn)

w = make_wcs(nhidata_fn)

RAstart = 15 # degrees, in Lazarian paper
RAstop = 35
DECstart = 4
DECstop = 16
x_start, y_start = radec_to_xy(RAstop, DECstart, w)
x_stop, y_stop = radec_to_xy(RAstart, DECstop, w)

print(x_start, y_start, x_stop, y_stop)

xcut_nhi_hdr, xcut_nhi_data = xcutout_data(nhi_data, nhi_hdr, xstart=x_start, xstop=x_stop)
ycut_nhi_hdr, ycut_nhi_data = ycutout_data(xcut_nhi_data, xcut_nhi_hdr, ystart=y_start, ystop=y_stop)


Qdata_fn = "/Volumes/DataDavy/GALFA/DR2/FullSkyRHT/new_thetarht_maps/QRHT_coadd_974_1069.fits"
Udata_fn = "/Volumes/DataDavy/GALFA/DR2/FullSkyRHT/new_thetarht_maps/URHT_coadd_974_1069.fits"
Qdata = fits.getdata(Qdata_fn)
Qhdr = fits.getheader(Qdata_fn)
Udata = fits.getdata(Udata_fn)
Uhdr = fits.getheader(Udata_fn)

xcut_U_hdr, xcut_U_data = xcutout_data(Udata, Uhdr, xstart=x_start, xstop=x_stop)
ycut_U_hdr, ycut_U_data = ycutout_data(xcut_U_data, xcut_U_hdr, ystart=y_start, ystop=y_stop)
xcut_Q_hdr, xcut_Q_data = xcutout_data(Qdata, Qhdr, xstart=x_start, xstop=x_stop)
ycut_Q_hdr, ycut_Q_data = ycutout_data(xcut_Q_data, xcut_Q_hdr, ystart=y_start, ystop=y_stop)

thetamap = 0.5*np.arctan2(ycut_U_data, ycut_Q_data)

ysize, xsize = thetamap.shape
X, Y = np.mgrid[0:xsize, 0:ysize]

skipint = 100

plt.figure()
plt.imshow(ycut_nhi_data, cmap="gray_r")

norm = np.sqrt(ycut_U_data**2+ycut_Q_data**2)
plotU = ycut_U_data/norm
plotQ = ycut_Q_data/norm
plt.quiver(X[::skipint, ::skipint], Y[::skipint, ::skipint], 
           plotU[::skipint, ::skipint], plotQ[::skipint, ::skipint], headaxislength=0, headlength=0)


#fig, ax = plt.subplots()
#im = ax.imshow(p, extent=[x.min(), x.max(), y.min(), y.max()])
#ax.quiver(x[skip], y[skip], dx[skip], dy[skip])


    