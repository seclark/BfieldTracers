from __future__ import division
import numpy as np
import cPickle
import matplotlib
matplotlib.rcParams['backend'] = "Qt4Agg"
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import wcs
import copy
from scipy import ndimage

import sys 
sys.path.insert(0, '../ACTPol/code')
import foreground_tools

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
    
def smooth_overnans(map, sig = 15):

    """
    Takes map with nans, etc set to 0
    """
    
    mask = np.ones(map.shape, np.float_)
    mask[np.isnan(map)] = 0
    
    map_zeroed = copy.copy(map)
    map_zeroed[mask == 0] = 0
    
    blurred_map = ndimage.gaussian_filter(map_zeroed, sigma=sig)
    blurred_mask = ndimage.gaussian_filter(mask, sigma=sig)
    
    map = blurred_map / blurred_mask
  
    return map
    
def make_SC_241_fits_QU():
    QRHT, URHT, PRHT, theta_rht, int_rhtunsmoothed, QRHTsq, URHTsq = foreground_tools.get_QU_RHT_corrected(region = "SC_241", wlen = 75, smr = 15, smoothRHT = False, sigma = 0, QUmean = False, bwrm = True, galfapixcorr = True, intRHTcorr = False)

    hdr = fits.getheader('/Volumes/DataDavy/GALFA/SC_241/cleaned/SC_241.66_28.675.best.fits')
    
    fits.writeto("/Volumes/DataDavy/GALFA/SC_241/thetarht_maps/QRHT_SC_241.66_28.675.best_w75_s15_t70_smoothRHT_False_bwrm_galfapixcorr.fits", QRHT, header=hdr)
    fits.writeto("/Volumes/DataDavy/GALFA/SC_241/thetarht_maps/URHT_SC_241.66_28.675.best_w75_s15_t70_smoothRHT_False_bwrm_galfapixcorr.fits", URHT, header=hdr)
    
    
nhidata_fn = "/Volumes/DataDavy/GALFA/DR2/NHIMaps/GALFA-HI_VLSR-036+0037kms_NHImap_noTcut.fits"
nhi_hdr = fits.getheader(nhidata_fn)
nhi_data = fits.getdata(nhidata_fn)

w = make_wcs(nhidata_fn)

RAstart = 215#15 # degrees, in Lazarian paper
RAstop = 260#35
DECstart = 20#4
DECstop = 36#16
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

planckQUdata_fn = "/Volumes/DataDavy/Planck/HFI_SkyMap_353_2048_R2.00_full_Equ_QUprojected_GALFAallsky.fits"
planckQU_data = fits.getdata(planckQUdata_fn)
planckQU_hdr = fits.getheader(planckQUdata_fn)
planckQdata = -planckQU_data[0, :, :].T
planckUdata = -planckQU_data[1, :, :].T
xcut_Q_hdr, xcut_Q_data = xcutout_data(planckQdata, planckQU_hdr, xstart=x_start, xstop=x_stop)
ycut_Q_hdr, PlanckQ = ycutout_data(xcut_Q_data, xcut_Q_hdr, ystart=y_start, ystop=y_stop)
xcut_U_hdr, xcut_U_data = xcutout_data(planckUdata, planckQU_hdr, xstart=x_start, xstop=x_stop)
ycut_U_hdr, PlanckU = ycutout_data(xcut_U_data, xcut_U_hdr, ystart=y_start, ystop=y_stop)

smoothsig = 120
smoothQ = smooth_overnans(ycut_Q_data, smoothsig)
smoothU = smooth_overnans(ycut_U_data, smoothsig)
smoothPlanckQ = smooth_overnans(PlanckQ, smoothsig)
smoothPlanckU = smooth_overnans(PlanckU, smoothsig)

thetamap = 0.5*np.arctan2(smoothU, smoothQ)

ysize, xsize = thetamap.shape
Y, X = np.mgrid[0:ysize, 0:xsize]

skipint = 500

plt.figure()
plt.imshow(ycut_nhi_data, cmap="gray_r")

pivot='mid'

norm = np.sqrt(smoothU**2+smoothQ**2)
plotU = smoothU/norm
plotQ = smoothQ/norm
#plt.quiver(X[::skipint, ::skipint], Y[::skipint, ::skipint], 
#           2*plotQ[::skipint, ::skipint], 2*plotU[::skipint, ::skipint], headaxislength=0, headlength=0, pivot=pivot)

#plt.quiver(X[::skipint, ::skipint], Y[::skipint, ::skipint], angles=np.degrees(thetamap[::skipint, ::skipint]), 
#           headaxislength=0, headlength=0)

plancknorm = np.sqrt(smoothPlanckU**2+smoothPlanckQ**2)
plotPlanckU = smoothPlanckU/plancknorm
plotPlanckQ = smoothPlanckQ/plancknorm
thetaplanck = 0.5*np.arctan2(smoothPlanckU, smoothPlanckQ)
#plt.quiver(X[::skipint, ::skipint], Y[::skipint, ::skipint], 
#           plotPlanckQ[::skipint, ::skipint], plotPlanckU[::skipint, ::skipint], headaxislength=0, headlength=0, pivot=pivot, color="yellow")
#X = X.ravel()
#Y = Y.ravel()
#plotPlanckQ = plotPlanckQ.ravel()
#plotPlanckU = plotPlanckU.ravel()
#plt.quiver(X[::skipint], Y[::skipint], 
#           plotPlanckQ[::skipint], plotPlanckU[::skipint], headaxislength=0, headlength=0, pivot=pivot, color="yellow")
plt.quiver(X[::skipint, ::skipint], Y[::skipint, ::skipint], X[::skipint, ::skipint], Y[::skipint, ::skipint], 
           thetaplanck[::skipint, ::skipint], headaxislength=0, headlength=0, pivot=pivot, color="yellow")



#fig, ax = plt.subplots()
#im = ax.imshow(p, extent=[x.min(), x.max(), y.min(), y.max()])
#ax.quiver(x[skip], y[skip], dx[skip], dy[skip])


    