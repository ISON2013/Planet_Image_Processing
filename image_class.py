#!/usr/bin/env python
# coding: utf-8

# In[1]:


from func.MapFunc import *
from tools.mapping import *
from tools.image import *
from tools.file import *
from astropy.time import Time
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
#get_ipython().run_line_magic('matplotlib', 'inline')
#plt.ion()

class image():
    
    def __init__(self, filepath):
        self.filepath = filepath
        self.img, self.hdr, self.scihdr = open_img(self.filepath) #opens the image
        
        if 'drz' in filepath:

            self.instrument = self.hdr['INSTRUME'] #pulls important values from fits header
            self.expstart = self.hdr['EXPSTART']
            self.exptime = self.hdr['EXPTIME']
            if 'WFC3' in self.instrument:
                self.filter = self.hdr['FILTER']
            else: 
                self.filter = self.hdr['FILTNAM1']

            self.orientat = self.scihdr['ORIENTAT'] #pulls important values from sci header
            self.identifier = self.scihdr['EXPNAME']
            self.photmode = self.scihdr['PHOTMODE']
            self.photflam = self.scihdr['PHOTFLAM']
            self.photplam = self.scihdr['PHOTPLAM']
            self.photbw = self.scihdr['PHOTBW']
            self.photzpt = self.scihdr['PHOTZPT']

            self.eph = get_ephemerides(target = 'neptune', loc='hst', epoch=self.expstart) #ephemeris from JPL Horizons

            self.nepdist = float(self.eph['r']) #pulls important values from ephemeris
            self.ang_width = float(self.eph['ang_width'])
            self.npole_ang = float(self.eph['NPole_ang'])
            self.obslong = float(self.eph['PDObsLon'])
            self.obslat = float(self.eph['PDObsLat'])
            self.sol_ang = -float(self.eph['alpha'])

            self.raw = self.img #saves a raw version of the image to rollback if needed.
        
        if 'WFC3' in self.instrument: #sets parameters which rely on other parameters
            self.platescale = 0.04
        elif 'WFPC2' in self.instrument:
            self.platescale = 0.05
        else:
            self.platescale = 0.04   #assumes an image without a valid instrument flag is a WFC3 image
        self.IF_Factor = 1
        self.planicentre = [0,0,0]
        self.radius = self.ang_width/(2*self.platescale)
        self.angcorr = -1*(float(self.orientat)-self.npole_ang)+180
        while self.angcorr >= 360 or self.angcorr < 0:
            if self.angcorr >=360:
                self.angcorr = self.angcorr - 360
            elif self.angcorr < 0:
                self.angcorr = self.angcorr + 360
        self.y_wide_crop = self.img.shape[0]
        self.x_wide_crop = self.img.shape[1]
        self.img_time = Time(float(self.expstart),format='mjd',scale='utc')
            
    def show(self,cmap = 'viridis', cmin=None, cmax=None, title=None, gridded=False):
        if cmin == None: #validates inputs
            cmin = np.nanmin(self.img)
        else:
            cmin = float(cmin)
        if cmax == None:
            cmax = np.nanmax(self.img)
        else:
            cmax = float(cmax)
        
        if title == None:
            title = self.identifier
        else:
            title = str(title)
        
        if not gridded:
            plt.figure(figsize=(5,5))
        plt.imshow(self.img)
        plt.colorbar()
        plt.clim(cmin,cmax)
        plt.title(title)
        if not gridded:
            plt.show()
            
    def trim_edge(self,border = -5):
        if border == -5:
            border = int(self.img.shape[0]*0.06)
        self.img = self.img[border:(self.img.shape[0]-border),border:(self.img.shape[1]-border)]
        self.trimmed = self.img
    
    def planicent(self, croppix = 400,Debug=False, accuracy = 10, returnval = False):
        print(self.identifier)
        if np.max(self.img) < 1:
            Threshold = 0.01
        else:
            Threshold = 5
        
        cropimg = self.img
        corner_trail_0 = []
        corner_trail_1 = []
        
        croploop = False
        imgcent = centroid_com(cropimg)
        
        croppix0 = 1e7
        croppix1 = 1e7
        if (imgcent[0]-croppix)<0:
            croppix0 = int(imgcent[0])
        elif (imgcent[0]+croppix)>cropimg.shape[0]:
            croppix0 = int(cropimg.shape[0]-imgcent[0])
        if (imgcent[1]-croppix)<0:
            croppix1 = int(imgcent[1])
        elif (imgcent[1]+croppix)>cropimg.shape[1]:
            croppix1 = int(cropimg.shape[1]-imgcent[1])
        if (croppix0 or croppix1) < croppix:
            if croppix0 < croppix1:
                croppix = croppix0
            else:
                croppix=croppix1
        else:
            pass
        
        while croploop == False:
            if Debug:
                plt.figure(figsize=(5,5))
                plt.imshow(cropimg)
                plt.axhline(imgcent[1])
                plt.axvline(imgcent[0])
                plt.colorbar()
                plt.show()
                
            if cropimg[int(imgcent[1]),int(imgcent[0])] < Threshold:
                corner_trail_0.append(imgcent[1]-croppix)
                corner_trail_1.append(imgcent[0]-croppix)
                cropimg = cropimg[(int(imgcent[1]-croppix)):(int(imgcent[1]+croppix)),(int(imgcent[0]-croppix)):(int(imgcent[0]+croppix))]
                croppix0 = 1e7
                croppix1 = 1e7
                if (imgcent[0]-croppix*.9)<0:
                    croppix0 = int(imgcent[0])
                elif (imgcent[0]+croppix*.9)>cropimg.shape[0]:
                    croppix0 = int(cropimg.shape[0]-imgcent[0])
                if (imgcent[1]-croppix*.9)<0:
                    croppix1 = int(imgcent[1])
                elif (imgcent[1]+croppix*.9)>cropimg.shape[1]:
                    croppix1 = int(cropimg.shape[1]-imgcent[1])
                if (croppix0 or croppix1) < croppix*.9:
                    if croppix0 < croppix1:
                        croppix = croppix0
                    else:
                        croppix=croppix1
                else:
                    croppix = croppix*.9
                imgcent = centroid_com(cropimg)
            else:
                croploop = True
        
        croppix0 = 1e7
        croppix1 = 1e7
        if (imgcent[0]-3*self.radius)<0:
            croppix0 = int(imgcent[0])
        elif (imgcent[0]+3*self.radius)>cropimg.shape[0]:
            croppix0 = int(cropimg.shape[0]-imgcent[0])
        if (imgcent[1]-3*self.radius)<0:
            croppix1 = int(imgcent[1])
        elif (imgcent[1]+3*self.radius)>cropimg.shape[1]:
            croppix1 = int(cropimg.shape[1]-imgcent[1])
        
        if (croppix0 or croppix1) < int(3*self.radius):
            if croppix0 < croppix1:
                croppix = croppix0
            else:
                croppix=croppix1
        else:
            croppix = int(3*self.radius)
            
        corner_trail_0.append(imgcent[1]-croppix)
        corner_trail_1.append(imgcent[0]-croppix)
        cropimg = cropimg[(int(imgcent[1]-croppix)):(int(imgcent[1]+croppix)),(int(imgcent[0]-croppix)):(int(imgcent[0]+croppix))]
        imgcent = get_disc(cropimg,r0=self.radius)
        
        if Debug:
            plt.figure(figsize=(5,5))
            plt.imshow(cropimg)
            plt.axhline(imgcent[1])
            plt.axvline(imgcent[0])
            plt.show()
        
        centrey = sum(corner_trail_0)+imgcent[1]
        centrex = sum(corner_trail_1)+imgcent[0]
        planicentre = [centrey,centrex,self.radius]
        
        #once confirmed planicentre on disc, finds limb of planet on a Factor 5 interpolated image
        img_interp, plani_interp = interp(self.img, planicentre, accuracy)
        yline = []
        for i in range(img_interp.shape[0]):
            yline.append(img_interp[i,int(plani_interp[1])])
        ygrad = np.gradient(yline)
        ylimb = [0,0]
        ylimb[0] = np.mean(np.where(ygrad == np.max(ygrad))[0])
        ylimb[1] = np.mean(np.where(ygrad == np.min(ygrad))[0])
        
        xline = []
        for i in range(img_interp.shape[1]):
            xline.append(img_interp[int(plani_interp[0]),i])
        xgrad = np.gradient(xline)
        xlimb = [0,0]
        xlimb[0] = np.mean(np.where(xgrad == np.max(xgrad))[0])
        xlimb[1] = np.mean(np.where(xgrad == np.min(xgrad))[0])
        
        if Debug:
            plt.figure(figsize=(15,3))
            plt.title('Y-axis')
            plt.axvline(ylimb[0], color='red')
            plt.axvline(ylimb[1], color='limegreen')
            plt.plot(np.linspace(0, len(ygrad), len(ygrad)), ygrad*4)
            plt.plot(np.linspace(0, len(yline), len(yline)), yline)
            plt.show()
            
            plt.figure(figsize=(15,3))
            plt.title('X-axis')
            plt.axvline(xlimb[0], color='red')
            plt.axvline(xlimb[1], color='limegreen')
            plt.plot(np.linspace(0, len(xgrad), len(xgrad)), xgrad*4)
            plt.plot(np.linspace(0, len(xline), len(xline)), xline)
            plt.show()
        
        plani_interp = [np.mean(ylimb),np.mean(xlimb),plani_interp[2]]
        if Debug:
            img_interp, y, x = circ(img_interp,plani_interp[2])
            plt.figure(figsize=(10,10))
            plt.imshow(img_interp)
            plt.axhline(plani_interp[0],color='limegreen')
            plt.axvline(plani_interp[1],color='limegreen')
            plt.plot(x+plani_interp[1], y+plani_interp[0],color='limegreen')
            plt.title('Interpolated, limb fitted planicentre')
            plt.show()

        self.planicentre = [plani_interp[0]/accuracy,plani_interp[1]/accuracy,self.radius]
        if Debug:
            plt.figure(figsize=(5,5))
            plt.imshow(self.img)
            plt.axhline(self.planicentre[0])
            plt.axvline(self.planicentre[1])
            plt.title('Final Planicentre')
            plt.show()

        print(self.identifier+': Planicentre Determined to be: '+str(self.planicentre))

        img_cent = [self.img.shape[0]/2,self.img.shape[1]/2]
        diff0 = self.planicentre[0]-img_cent[0]
        diff1 = self.planicentre[1]-img_cent[1]
        self.y_wide_crop = img_cent[0]-abs(diff0)
        self.x_wide_crop = img_cent[1]-abs(diff1)
        if returnval:
            return [self.planicentre[0], self.planicentre[1]]
        else:
            pass
        
    def back_sub(self, corner = 'topleft', bbox_factor = 0.3, Debug = False):
        size_y = int(self.img.shape[0]*bbox_factor)
        size_x = int(self.img.shape[1]*bbox_factor)
        
        if 'top' in corner:
            ybounds = [0,size_y]
        elif 'bottom' in corner:
            ybounds = [self.img.shape[0]-size_y,self.img.shape[0]]
        else:
            raise ValueError('Invalid entry string, must contain either "top" or "bottom".')
            
        if 'left' in corner:
            xbounds = [0,size_x]
        elif 'right' in corner:
            xbounds = [self.img.shape[1]-size_x,self.img.shape[1]]
        else:
            raise ValueError('Invalid entry string, must contain either "left" or "right".')
        
        if Debug:
            plt.figure(figsize=(10,5))
            plt.imshow(self.img)
            plt.axhline(ybounds[0], color = 'green')
            plt.axhline(ybounds[1], color = 'green')
            plt.axvline(xbounds[0], color = 'green')
            plt.axvline(xbounds[1], color = 'green')
            plt.show()
        
        background = np.mean(self.img[ybounds[0]:ybounds[1], xbounds[0]:xbounds[1]])
        print('Background brightness of '+str(background)+' subtracted.')
        self.img = self.img-background
        
    def calcIF(self, flux = None):
        if flux == None:
            lookup = self.instrument.lower()+self.filter.lower()
            
            with open('C:/Users/jamie/Documents/Planet_Image_Processing/ancilliary/filter_fluxes/filter_solar_fluxes.txt', 'r') as fileobj:
                for line in fileobj:
                    if lookup in line:
                        flux = float(str(line).split(',')[2])

        self.IF, badscale = calc_IF(self.nepdist, self.photflam, flux, PlateScale = self.platescale)
    
    def calcPhys(self):
        pass
    
    def applyFactor(self, factor):
        self.img = self.img*factor
        
    def calib(self, units='IF', back_loc = 'topleft', flux = None):
        self.back_sub(back_loc)
        if units == 'IF':
            self.calcIF(flux)
            self.applyFactor(self.IF)
        elif units == 'phys':
            self.calcPhys()
            self.applyFactor(factor)
    
    def shift(self, loc = None, accuracy = 10):
    
        #interpolate image to get to the accuracy needed
        img_interp, plani_interp = interp(self.img, self.planicentre, accuracy) 
        if loc == None:
            loc = [img_interp.shape[0]/2, img_interp.shape[1]/2]
        else:
            loc = [loc[0]*accuracy, loc[1]*accuracy]
        
        centre = [plani_interp[0],plani_interp[1]] #perform the centering
        diff = np.array(loc)-np.array(centre)
        diff = diff.astype('int32')
        diff = diff.tolist()
        img_interp = np.roll(img_interp, diff, axis=[0,1])
        plani_interp = [loc[0],loc[1],plani_interp[2]]
        
        #'deinterpolate' the image back to prior resolution
        self.img, self.planicentre = interp(img_interp, plani_interp, 1/accuracy)
        
    def crop(self, crop = 'wide'):
        if 'wide' in crop:
            self.img = self.img[round(self.planicentre[0]-self.y_wide_crop):round(self.planicentre[0]+self.y_wide_crop),round(self.planicentre[1]-self.x_wide_crop):round(self.planicentre[1]+self.x_wide_crop)]
        elif 'superzoom' in crop:
            self.img = self.img[round(self.planicentre[0]-self.radius*1.3):round(self.planicentre[0]+self.radius*1.3),round(self.planicentre[1]-self.radius*1.3):round(self.planicentre[1]+self.radius*1.3)]
        elif 'close' in crop:
            self.img = self.img[round(self.planicentre[0]-self.radius*2):round(self.planicentre[0]+self.radius*2),round(self.planicentre[1]-self.radius*2):round(self.planicentre[1]+self.radius*2)]
        else:
            raise ValueError('Invalid entry string, must contain either "wide", "close" or "superzoom".')
        self.planicentre = [self.img.shape[0]/2, self.img.shape[1]/2, self.planicentre[2]]
    
    def rotate(self):
        self.img = ndimage.rotate(self.img, self.angcorr, reshape=False)
    
    def minnaert(self, k, img, ex_ang_obs=70, ex_ang_sun=70):
        k = float(k)
        thetainc = self.sol_ang
        thetainc %= 360

        #find the theta and phi coordinate images for the disc
        x, y, z = img_to_xyz(img.shape, self.planicentre[0], self.planicentre[1], self.radius)
        theta, phi = xyz_to_thetaphi(x, y, z)
        theta_i, phi_i = xyz_to_thetaphi(x,y,z,theta=thetainc)

        #get the mu parameter image
        mu = np.cos(theta*(np.pi/180))
        mu0 = np.cos(theta_i*(np.pi/180))

        #apply the correction to the image
        img = img/((mu**(k-1))*(mu0**k))
        img[theta>ex_ang_obs] = np.nan #trims the extreme edges of the planetary disc off
        img[theta_i>ex_ang_sun] = np.nan #trims the night side of the planet off
        
        return img
    
    def change_radius(self, radius_goal = 200, images = ['img','background','cloud']):
        factor = radius_goal/self.radius
        if 'img' in images:
            self.img, planicentre_img = interp(self.img, self.planicentre, factor)
        if 'background' in images:
            self.background, planicentre = interp(self.background, self.planicentre, factor)
        if 'cloud' in images:
            self.cloud_fill, planicentre = interp(self.cloud_fill, self.planicentre, factor)
        
        self.planicentre = planicentre_img
        self.radius = radius_goal
    
    def show_map(self, map_type = 'img', cmap = 'viridis', lvls=100, cmin=None, cmax=None, title=None, gridded=False):
        if cmin == None: #validates inputs
            if map_type == 'img':
                cmin = np.nanmin(self.mapped_img[0])
            elif map_type == 'background':
                cmin = np.nanmin(self.mapped_back[0])
            elif map_type == 'cloud':
                cmin = np.nanmin(self.mapped_cloud[0])
        else:
            cmin = float(cmin)
        if cmax == None:
            if map_type == 'img':
                cmax = np.nanmax(self.mapped_img[0])
            elif map_type == 'background':
                cmax = np.nanmax(self.mapped_back[0])
            elif map_type == 'cloud':
                cmax = np.nanmax(self.mapped_cloud[0])
        else:
            cmax = float(cmax)
        
        if title != None:
            title = str(title)
        elif map_type == 'img':
            title = str(self.identifier)+': Image Mapped'
        elif map_type == 'background':
            title = str(self.identifier)+': Background Image Mapped'
        elif map_type == 'cloud':
            title = str(self.identifier)+': Cloud Image Mapped'
        
        if not gridded:
            plt.figure(figsize=(12.5,5))
        if map_type == 'img':
            plt.contourf(self.mapped_img[1], self.mapped_img[2], self.mapped_img[0],lvls,cmap=cmap)
        elif map_type == 'background':
            plt.contourf(self.mapped_back[1], self.mapped_back[2], self.mapped_back[0],lvls,cmap=cmap)
        elif map_type == 'cloud':
            plt.contourf(self.mapped_cloud[1], self.mapped_cloud[2], self.mapped_cloud[0],lvls,cmap=cmap)
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.title(title)
        plt.colorbar()
        plt.clim(cmin,cmax)
        if not gridded:
            plt.show()
    
    def map_image(self, maps = ['img','background','cloud'], showMaps = True):
        # create maps of images
        if 'img' in maps:
            self.mapped_img = map_observation(self.img, self.img.shape[0]/2, self.img.shape[1]/2, self.radius, 180, self.obslong, self.obslat)
        if 'background' in maps:
            self.mapped_back = map_observation(self.background, self.background.shape[0]/2, self.background.shape[1]/2, self.radius, 180, self.obslong, self.obslat)
        if 'cloud' in maps:
            self.mapped_cloud = map_observation(self.cloud_fill, self.cloud_fill.shape[0]/2, self.cloud_fill.shape[1]/2, self.radius, 180, self.obslong, self.obslat)
        
        if showMaps:
            plt_xdim = 12.5*len(maps)
            plt.figure(figsize=(plt_xdim,5))
            
            for i in range(len(maps)):
                ax = plt.subplot(1,len(maps),i+1)
                self.show_map(map_type = maps[i], gridded = True)
            plt.show()
    
    def split_img(self, Debug = False): #splits an image into 2 vertical halves, since planet is centred, also splits planet
        if self.img.shape[1]%2 == 0:
            deinterp = False
        else:
            self.img, self.planicentre = interp(self.img, self.planicentre, factor=2)
            deinterp = True
        
        left_img = self.img[:,0:int(self.img.shape[1]/2)]
        right_img = self.img[:,int(self.img.shape[1]/2):int(self.img.shape[1])]
        
        if deinterp:
            self.img, self.planicentre = interp(self.img, self.planicentre, factor = 0.5)
            left_img, ignore = interp(left_img, self.planicentre, factor=0.5)
            right_img, ignore = interp(right_img, self.planicentre, factor=0.5)
        if Debug:
            self.colorbarlim = [np.min(self.img), np.max(self.img)]
            plt.figure(figsize=(10,10))
            ax_left = plt.subplot(121)
            plt.imshow(left_img)
            plt.colorbar()
            plt.clim(self.colorbarlim[0],self.colorbarlim[1])
            ax_right = plt.subplot(122)
            plt.imshow(right_img)
            plt.colorbar()
            plt.clim(self.colorbarlim[0],self.colorbarlim[1])
            plt.show()
        return left_img, right_img
    
print('image class initialised.')
