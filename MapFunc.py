# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 14:46:03 2021

@author: Jamie MacQuillin
"""
#Import all the necessary libraries
import matplotlib.pyplot as plt
import numpy as np
from photutils import centroid_com
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.time import Time

from tools import *
from tools.mapping import *
from tools.image import *

"""

"""
def open_img(path, ext = 'sci'):
    #shield against non-string inputs
    path = str(path)
    ext_dict = {
        'sci': 1, 
        }
    if ext in ext_dict:
        ext = ext_dict[ext]
        
    
    #filename = get_pkg_data_filename(path)
    fitsdata = fits.open(path)
    imgdata = fits.getdata(path, ext=1)
    fitshdr = fitsdata[0].header
    scihdr = fitsdata[1].header
    return imgdata, fitshdr, scihdr   

"""
Functions which aid with the calibration of the images

- AU_to_km: Takes a distance in AU, converts it to km
                n : distance (km)

- calc_IF: Calculates an I/F factor for a non byte-scaled image, defaults
             based on Hubble WFC3/UVIS instrument
                r : distance between sun and planet (km)
         PHOTFLAM : Inverse sensitivity of the instrument erg/cm^-2/ang/electron
        SolarFlux : Flux density of the sun in that wavelength at Earth
        PlateScale: Plate Scale of the detector (arcsec)

- minnaert: Applies a Minnaert Correction to an image of a disc
"""
def AU_to_km(n):
    Dist = n*149598000
    return Dist

def calc_IF(r, PHOTFLAM, SolarFlux, PlateScale = 0.04, Filter = None):
    # Shielding from arrays of length 1
    dist = float(r)
    PHOTFLAM = float(PHOTFLAM)
    SolFlux = float(SolarFlux)
    PlateScale = float(PlateScale)
    
    if Filter == '467':
        #determining I/F Factor
        omega =((0.0495/3600)*(np.pi/180))**2 #3.76E-14
        fractop = ((PHOTFLAM*10)/omega)
        fracbottom = ((AU_to_km(1.0)/AU_to_km(dist))**2)*(SolFlux/np.pi)
        IF = fractop/fracbottom#PHOTFLAM*(10/omega)*(np.pi/SolFlux)*(dist**2)
    elif Filter == '547':
        #determining I/F Factor
        omega =((0.049/3600)*(np.pi/180))**2 #3.76E-14
        fractop = ((PHOTFLAM*10)/omega)
        fracbottom = ((AU_to_km(1.0)/AU_to_km(dist))**2)*(SolFlux/np.pi)
        IF = fractop/fracbottom#PHOTFLAM*(10/omega)*(np.pi/SolFlux)*(dist**2)
    elif Filter == '619':
        #determining I/F Factor
        omega =((0.0497/3600)*(np.pi/180))**2 #3.76E-14
        fractop = ((PHOTFLAM*10)/omega)
        fracbottom = ((AU_to_km(1.0)/AU_to_km(dist))**2)*(SolFlux/np.pi)
        IF = fractop/fracbottom#PHOTFLAM*(10/omega)*(np.pi/SolFlux)*(dist**2)
    else:
        #determining I/F Factor
        omega =((PlateScale/3600)*(np.pi/180))**2 #3.76E-14
        fractop = ((PHOTFLAM*10)/omega)
        fracbottom = ((AU_to_km(1.0)/AU_to_km(dist))**2)*(SolFlux/np.pi)
        IF = fractop/fracbottom#PHOTFLAM*(10/omega)*(np.pi/SolFlux)*(dist**2)
    return IF, PlateScale

def minnaert(img, cent, k, eph, ex_ang=70):
    """
    Function which applies a Minnaert Function Correction to an image, removing limb darkening
    Inputs are:
    - img: the cropped image of the planetary disc
    - cent: the planicentre and radius of the disc, as given by the get_disc routine from Oliver King
    - k: the minnaert parameter to be applied to the image
    - eph: the ephemeris object corresponding to the image
    - ex_ang: the emission angle from which to exclude the outer edges of the disc, default value 70 degrees
    
    Outputs are:
    - img_corr: the Minnaert Corrected image of the disc
    """
    #Shielding from unexpected forms of inputs for k
    k = float(k)
    thetainc = -float(eph['alpha'])
    thetainc %= 360
    
    #find the theta and phi coordinate images for the disc
    x, y, z = img_to_xyz(img.shape, cent[0], cent[1], cent[2])
    theta, phi = xyz_to_thetaphi(x, y, z)
    theta_i, phi_i = xyz_to_thetaphi(x,y,z,theta=thetainc)
    
    #get the mu parameter image
    mu = np.cos(theta*(np.pi/180))
    mu0 = np.cos(theta_i*(np.pi/180))
    
    #apply the correction to the image
    img_corr = img/((mu**(k-1))*(mu0**k))
    img_corr[theta>ex_ang] = np.nan #trims the extreme edges of the planetary disc off
    return img_corr

def vgr_minnaert(img, cent, k, eph, ex_ang_obs=70, ex_ang_sun=70):
    """
    Function which applies a Minnaert Function Correction to an image, removing limb darkening
    Inputs are:
    - img: the cropped image of the planetary disc
    - cent: the planicentre and radius of the disc, as given by the get_disc routine from Oliver King
    - k: the minnaert parameter to be applied to the image
    - eph: the ephemeris object corresponding to the image
    - ex_ang: the emission angle from which to exclude the outer edges of the disc, default value 70 degrees
    
    Outputs are:
    - img_corr: the Minnaert Corrected image of the disc
    """
    #Shielding from unexpected forms of inputs for k
    k = float(k)
    thetainc = -float(eph['alpha'])
    thetainc %= 360
    
    #find the theta and phi coordinate images for the disc
    x, y, z = img_to_xyz(img.shape, cent[0], cent[1], cent[2])
    theta, phi = xyz_to_thetaphi(x, y, z)
    theta_i, phi_i = xyz_to_thetaphi(x,y,z,theta=thetainc)
    
    #get the mu parameter image
    mu = np.cos(theta*(np.pi/180))
    mu0 = np.cos(theta_i*(np.pi/180))
    
    #apply the correction to the image
    img_corr = img/((mu**(k-1))*(mu0**k))
    img_corr[theta>ex_ang_obs] = np.nan #trims the extreme edges of the planetary disc off
    img_corr[theta_i>ex_ang_sun] = np.nan #trims the night side of the planet off
    return img_corr

"""
Functions which aid with navigation of the images

- circ: Cuts out a circle around the image
            img : the full image of the planet
            rad : the radius of the disc of the planet
            
- get_planicentre: Finds the planicentre of a given disc (starred inputs are required)
            img* : the full image of the planet
            rad* : the radius of the disc of the planet
            This routine can iterate if it does not find the planicentre after a default run through, to find the actual disk
            It does this by cropping the image around the first guess, and then running the guess process iteratively until it converges
            iterations : the number of times the iterative process should be run
            croppix : the starting number of pixels each side of the guess at the planicentre that should be left at iteration 1
            Debug : a boolean variable determining if the images are plotted at each iteration, used when setting iteration and croppix
        returns:
            DiscImg : cropped image of the disc (subarray of the imported img)
            planicentre2 : an object containing the 2 coordinates of the centre
                           and the radius of the disc
"""
def circ(img, rad):
    theta = np.linspace(0, 2*np.pi, 100)
    r = np.sqrt(rad**2)
    y = r*np.cos(theta)
    x = r*np.sin(theta)
    return(img, y, x)

def circ2(img, rad):
    theta = np.linspace(0, 2*np.pi, 100)
    r = np.sqrt(rad**2)
    y = r*np.cos(theta)
    x = r*np.sin(theta)
    return(img, x, y)

def get_planicentre(img, rad, iterations=1, croppix=400, Debug=False,endcrop=10,despike=True,crop=True):
    croppedimg = img
        
    for i in range(iterations):
        imgcent = centroid_com(croppedimg)
        
        if Debug:
            print(imgcent)
            plt.figure(figsize=(5,5))
            im5 = plt.imshow(croppedimg)
            plt.axhline(imgcent[0],color='red')
            plt.axvline(imgcent[1],color='red')
            plt.title('Partial Image '+str(i))
            plt.colorbar(im5)
            plt.show() 
        
        if crop:
            croppedimg = croppedimg[(int(imgcent[1]-croppix)):(int(imgcent[1]+croppix)),(int(imgcent[0]-croppix)):(int(imgcent[0]+croppix))]
            croppix = croppix/3
        if despike:
            croppedimg = image.exp_despike(croppedimg, nsigma=1)
       
    planicentre1 = centroid_com(croppedimg)
    #print('get_planicentre1'+str(planicentre1))
    if crop:
        DiscImg = croppedimg[(int(planicentre1[1])-(int(rad)+endcrop)):(int(planicentre1[1])+(int(rad)+endcrop)),(int(planicentre1[0])-(int(rad)+endcrop)):int(planicentre1[0])+((int(rad)+endcrop))]
    else:
        DiscImg = croppedimg
    planicentre2 = get_disc(DiscImg, r0=rad)
    #print('get_planicentre2:'+str(planicentre2))
    return DiscImg, planicentre2

def flipplanicentre(planicentre):
    planicentre0 = planicentre[0]
    planicentre1 = planicentre[1]
    return planicentre0, planicentre0

def get_vgr_planicentre(img, rad, pltscl, iterations=1, croppix=400, Debug=False):  
    planicentre = centroid_com(img)
    xline = []
    for i in range(len(img)):
        xline.append(img[int(planicentre[1]),i])
    xgrad = np.gradient(xline)
    xlimb = [0,0]
    xlimb[0] = int(np.where(xgrad == np.amax(xgrad))[0])
    xlimb[1] = int(np.where(xgrad == np.amin(xgrad))[0])
    
    if not 0.25<abs(xgrad[xlimb[1]])/abs(xgrad[xlimb[0]])<1.75:
        if xgrad[xlimb[1]]>xgrad[xlimb[0]]:
            xlimb[0] = xlimb[1]-2*rad
        else:
            xlimb[1] = xlimb[0]+2*rad
    
    yline = []
    for i in range(len(img[int(planicentre[0])])):
        yline.append(img[i,int(planicentre[1])])
    ygrad = np.gradient(yline)
    ylimb = [0,0]
    ylimb[0] = int(np.where(ygrad == np.amax(ygrad))[0])
    ylimb[1] = int(np.where(ygrad == np.amin(ygrad))[0])
    
    if not 0.25<abs(ygrad[ylimb[1]])/abs(ygrad[ylimb[0]])<1.75:
        if ygrad[ylimb[1]]>ygrad[ylimb[0]]:
            ylimb[0] = ylimb[1]-2*rad
        else:
            ylimb[1] = ylimb[0]+2*rad
    
    if Debug:
        plt.figure(figsize=(15,3))
        plt.title('X-axis')
        plt.axvline(xlimb[0], color='red')
        plt.axvline(xlimb[1], color='limegreen')
        plt.plot(np.linspace(0, len(xgrad), len(xgrad)), xgrad*4)
        plt.plot(np.linspace(0, len(xline), len(xline)), xline)
        plt.show()

        plt.figure(figsize=(15,3))
        plt.title('Y-axis')
        plt.axvline(ylimb[0], color='red')
        plt.axvline(ylimb[1], color='limegreen')
        plt.plot(np.linspace(0, len(ygrad), len(ygrad)), ygrad*4)
        plt.plot(np.linspace(0, len(yline), len(yline)), yline)
        plt.show()
    
    planicentre2 = [np.mean(xlimb),np.mean(ylimb),rad]
    
    img, y, x = circ(img,rad)#generates required variables to draw circle on the image
    
    if Debug:
        fig = plt.figure(figsize=(10,10))
        ax2 = plt.subplot(111)#122
        im2 = ax2.imshow(img, cmap='bone')
        plt.title('Planicentre Guesses')
        ax2.axhline(planicentre[1],color='red')
        ax2.axvline(planicentre[0],color='red')
        ax2.plot(x+planicentre[0], y+planicentre[1],color='red')
        ax2.axvline(xlimb[0],color='limegreen')
        ax2.axvline(xlimb[1],color='limegreen')
        ax2.axvline(np.mean(xlimb),color='limegreen')
        ax2.axhline(ylimb[0],color='limegreen')
        ax2.axhline(ylimb[1],color='limegreen')
        ax2.axhline(np.mean(ylimb),color='limegreen')
        ax2.plot(x+planicentre2[0], y+planicentre2[1],color='limegreen')
        plt.show()
    return img, planicentre2

def maptrim(mapped_img, Debug = False):
    longfact = (360/len(mapped_img[1]))
    latfact = (180/len(mapped_img[2]))
    
    #loop to find centroid of map
    trimlong = 90/longfact #variable which gives the total width of the mapped segment desired
    equator = int(90/latfact) #gives the equator location in map pixel coordinates
    
    limb = [0,0]
    equatorline = mapped_img[0][equator]
    equatorline[np.where(np.isnan(equatorline))] = 0
    equ_grad = np.gradient(equatorline)
    limb[0] = int(np.where(equ_grad == np.amax(equ_grad))[0])
    limb[1] = int(np.where(equ_grad == np.amin(equ_grad))[0])
    
    if Debug: 
        plt.figure(figsize=(15,7))
        plt.plot(mapped_img[1],equatorline)
        plt.plot(mapped_img[1],equ_grad)
        plt.axvline(limb[0]*longfact)
        plt.axvline(limb[1]*longfact)
        plt.show()

    centlong = (limb[0]+limb[1])/2
    if Debug:
        print(centlong)
    
    trimmed_img = [0,0,0]
    trimmed_img[0] = mapped_img[0][:,int(centlong)-int(trimlong/2):int(centlong)+int(trimlong/2)]
    trimmed_img[1] = mapped_img[1][int(centlong)-int(trimlong/2):int(centlong)+int(trimlong/2)]
    trimmed_img[2] = mapped_img[2]
    return trimmed_img

def readvgrlbl(productid):
    file = open('../ALL_GEOMED_IMAGES/'+productid+'_GEOMED.lbl', 'r')
    labelfile = file.readlines()
    file.close()
    
    ids = []
    for i in range(len(labelfile)):
        ids.append(str(labelfile[i][0:32]))
    
    IF_Loc = ids.index('  REFLECTANCE_SCALING_FACTOR    ')
    Epoch_Loc = ids.index('IMAGE_TIME                      ')
    pix_fov_Loc = ids.index('  HORIZONTAL_PIXEL_FOV          ')
    
    IF = float(labelfile[IF_Loc][34:])
    datetime = str(labelfile[Epoch_Loc][34:-1])
    pix_fov = float(labelfile[pix_fov_Loc][34:44])
    
    t = Time(datetime)
    epoch = t.mjd
    return IF, epoch, pix_fov

def interp(image, planicentre, factor):
    factor = float(factor)
    img_interp = interp_image(image, factor)
    planicentre_interp = [planicentre[0]*factor,planicentre[1]*factor,planicentre[2]*factor]
    return img_interp, planicentre_interp

pltscl_WFC3 = 0.04
pltscl_WFPC2 = 0.046