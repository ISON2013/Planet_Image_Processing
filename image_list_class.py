import fnmatch
import os
from tools.file import *
from tools.image import *
import astropy.io.fits as ast
from astropy.time import Time
import math
from func.image_class import *

def getListOfFiles(dirName):
    # create a list of file and sub directories 
    # names in the given directory 
    listOfFile = os.listdir(dirName)
    allFiles = list()
    # Iterate over all the entries
    for entry in listOfFile:
        # Create full path
        fullPath = os.path.join(dirName, entry)
        # If entry is a directory then get the list of files in this directory 
        if os.path.isdir(fullPath):
            allFiles = allFiles + getListOfFiles(fullPath)
        else:
            allFiles.append(fullPath)
    return allFiles

def filesearch(filestring):
    """Function which returns list of filepaths found by wildcard searching
    the filepaths given"""
    dirs = filestring.split('/')
    basestr = ''
    i=0
    while '*' not in dirs[i]:
        basestr = basestr+dirs[i]+'/'
        i += 1
    filelist = getListOfFiles(basestr)
    filelist = fnmatch.filter(filelist, filestring)
    return filelist

back_locations = ['topleft','topright','bottomleft','bottomright']

class image_list():
    
    def __init__(self):
        self.filepaths = []
        self.img_list = []
        self.img_ids = []
        self.avg_background = []
        
    def add_images(self, filepaths):
        self.filepaths = filepaths
        for elem in self.filepaths:
            self.img_list.append(image(elem))
        for elem in self.img_list:
            self.img_ids.append(elem.identifier)
    
    def show_grid(self,img_type = 'img',cmin = 0, cmax = None, cmap = 'viridis'):
        if cmax == None:
            cmax = 0
            for i in range(len(self.img_list)):
                if img_type == 'img':
                    tempmax = np.nanmax(self.img_list[i].img)
                elif img_type == 'background':
                    tempmax = np.nanmax(self.img_list[i].background)
                elif img_type == 'cloud':
                    tempmax = np.nanmax(self.img_list[i].cloud_fill)
                else:
                    raise ValueError('img_type must be "img". "background" or "cloud".')
                if tempmax > cmax:
                    cmax = tempmax
        
        if len(self.img_list)%5 == 0:
            rows = len(self.img_list)//5
        else:
            rows = (len(self.img_list)//5)+1
        plt_y_len = 5*(rows)
        
        plt.figure(figsize=(25,plt_y_len))
        for i in range(len(self.img_list)):
            ax = plt.subplot(rows,5,i+1)
            if img_type == 'img':
                    plt.imshow(self.img_list[i].img,cmap = cmap)
            elif img_type == 'background':
                plt.imshow(self.img_list[i].background,cmap = cmap)
            elif img_type == 'cloud':
                plt.imshow(self.img_list[i].cloud_fill,cmap = cmap)
            else:
                raise ValueError('img_type must be "img". "background" or "cloud".')
            plt.title(self.img_list[i].identifier)
            plt.colorbar()
            plt.clim(cmin,cmax)
        plt.show()
    
    def img_goodness(self, badimgs = []):
        if badimgs == []:
            self.show_grid()
            print('List bad image IDs. Press enter between each id. Type "done" to move on.')
            bad_id = input(':> ')

            while 'done' not in bad_id:
                location = self.img_ids.index(bad_id)
                del self.img_list[location]
                del self.img_ids[location]
                bad_id = input(':> ')
        else:
            for elem in badimgs:
                if 'done' in elem:
                    pass
                else:
                    location = self.img_ids.index(elem)
                    del self.img_list[location]
                    del self.img_ids[location]
                
        print('Bad Images Removed.')
    
    def trim_edge_all(self, border = -5):
        for elem in self.img_list:
            elem.trim_edge(border)
    
    def planicent_all(self, croppix = 400, accuracy = 10, Debug = False):
        for elem in self.img_list:
            elem.planicent(croppix, Debug, accuracy)
    
    def calib_all(self, units='IF', back_loc = 'topleft', flux = None):
        for elem in self.img_list:
            elem.calib(units, back_loc, flux)
    
    def shift_all(self, loc=None, accuracy=10):
        for elem in self.img_list:
            elem.shift(loc, accuracy)
            if loc == None:
                print(str(elem.identifier)+' moved to centre.')
            else:
                print(str(elem.identifier)+' moved to '+print(str(loc)))
    
    def normalise_radius(self):
        max_radius = 0
        for elem in self.img_list:
            if elem.radius>max_radius:
                max_radius = elem.radius
        goal_rad = math.ceil(max_radius)
        
        for elem in self.img_list:
            elem.change_radius(radius_goal = goal_rad, images = ['img'])
            print(str(elem.identifier)+' radius changed to '+str(elem.radius))
    
    def crop_all(self,crop = 'wide'):
        for elem in self.img_list:
            elem.crop(crop=crop)
        if crop == 'wide':
            sizes0 = []
            sizes1 = []
            for elem in self.img_list:
                sizes0.append(elem.img.shape[0])
                sizes1.append(elem.img.shape[1])
            print(sizes0)
            print(sizes1)
    
    def rotate_all(self):
        for elem in self.img_list:
            elem.rotate()
    
    def isolate_backgrounds(self, showBackground = True, averageCutoff = 1, cloudAlbedo = 0.9, Debug=False):
        colorbarlim = [0,0]
        self.mirror_imgs = []
        
        for elem in self.img_list:
            tempmax = np.max(elem.img)
            if colorbarlim[1]<tempmax:
                colorbarlim[1] = tempmax
            templeft, tempright = elem.split_img()
            self.mirror_imgs.append(np.concatenate((templeft,np.fliplr(templeft)),axis=1))
            self.mirror_imgs.append(np.concatenate((np.fliplr(tempright),tempright),axis=1))
        
        max_y = 0
        max_x = 0
        for elem in self.mirror_imgs:
            y = elem.shape[0]
            x = elem.shape[1]
            if y>max_y:
                max_y = y
            if x>max_x:
                max_x = x
        #print('y = '+str(max_y))
        #print('x = '+str(max_x))
        #print(len(self.mirror_imgs))
        
        for k in range(len(self.mirror_imgs)):
            if (self.mirror_imgs[k].shape[0] != max_y) or (self.mirror_imgs[k].shape[1] != max_x):
                padarray = np.zeros([max_y, max_x])
                for i in range(self.mirror_imgs[k].shape[0]):
                    for j in range(self.mirror_imgs[k].shape[1]):
                        padarray[i,j] = self.mirror_imgs[k][i,j] 
                self.mirror_imgs[k] = padarray
        
        cutoff = int(len(self.mirror_imgs)*averageCutoff)
        self.avg_background = np.zeros([max_y, max_x])
        
        if Debug:
            #show grid of mirror images
            cmin = 0
            cmax = 0
            for i in range(len(self.mirror_imgs)):
                tempmax = np.max(self.mirror_imgs[i])
                if tempmax > cmax:
                    cmax = tempmax
            rows = (len(self.mirror_imgs)//5)+1
            plt_y_len = 5*(rows)
            plt.figure(figsize=(25,plt_y_len))
            for i in range(len(self.mirror_imgs[i])):
                ax = plt.subplot(rows,5,i+1)
                plt.imshow(self.mirror_imgs[i])
                plt.colorbar()
                plt.clim(cmin,cmax)
            plt.show()
        
        for y in range(max_y):
            for x in range(max_x):
                pixels = []
                for k in range(len(self.mirror_imgs)):
                    pixels.append(self.mirror_imgs[k][y,x])
                    #print(k)
                pixels.sort()
                pixels = pixels[0:cutoff]
                self.avg_background[y,x] = np.mean(pixels)
        
        if showBackground:
            plt.figure(figsize=(5,5))
            plt.imshow(self.avg_background)
            plt.title('Average Background Image')
            plt.colorbar()
            plt.show()
        
        for elem in self.img_list:
            print(elem.identifier+': Background and Cloud Separated')
            elem.cloud_fill = (elem.img-self.avg_background)/(cloudAlbedo-self.avg_background)
            elem.background = (elem.img-(cloudAlbedo*elem.cloud_fill))/(1-elem.cloud_fill)
    
    def change_radii_all(self, goalrad = 200, images = ['img','background','cloud']):
        for elem in self.img_list:
            elem.change_radius(goalrad, images)
        
    def minnaert_all(self, k, images = ['img','background','cloud']):
        for elem in self.img_list:
            if 'img' in images:
                elem.img = elem.minnaert(k, elem.img)
            if 'background' in images:
                elem.background = elem.minnaert(k, elem.background)
            if 'cloud' in images:
                elem.cloud_fill = elem.minnaert(k, elem.cloud_fill)
            print(str(elem.identifier)+' minnaert corrected, with minnaert factor k='+str(k))
    
    def map_all(self, images = ['img','background','cloud'], showMaps = True):
        for elem in self.img_list:
            elem.map_image(images, showMaps)
            print(str(elem.identifier)+' mapped successfully')
    
    def raw_to_png(self, k, wavlen, badimgs = [],pngmin = 0, pngmax = 1,averageCutoff = 0.5):
        self.trim_edge_all()
        self.img_goodness(badimgs = badimgs)
        self.planicent_all()
        self.calib_all()
        self.shift_all()
        self.normalise_radius()
        self.crop_all(crop='close')
        self.rotate_all()
        #self.planicent_all()
        self.write_img_to_png('C:/Users/jamie/Documents/Planet_Image_Processing/processed_data/pngs/'+self.img_list[0].filter+'/', wavlen = wavlen, k=0,cmin=pngmin, cmax = pngmax)
        self.isolate_backgrounds(showBackground = False, averageCutoff = averageCutoff)
        self.write_to_fits(filepath = 'C:/Users/jamie/Documents/Planet_Image_Processing/processed_data/fits/',crop_of_img='close')
        
    def raw_to_map(self, k, wavlen, badimgs = [], averageCutoff = 0.5, show_results = False,pngmin = 0, pngmax = 1):
        self.trim_edge_all()
        self.img_goodness(badimgs = badimgs)
        self.planicent_all()
        self.calib_all()
        self.shift_all()
        self.normalise_radius()
        self.crop_all(crop='superzoom')
        self.rotate_all()
#         self.planicent_all()
#         self.shift_all()
                
        if show_results:
            self.isolate_backgrounds(averageCutoff = averageCutoff)
        else:
            self.isolate_backgrounds(showBackground = False, averageCutoff = averageCutoff)
        self.write_to_fits(filepath = 'C:/Users/jamie/Documents/Planet_Image_Processing/processed_data/fits/',crop_of_img='superzoom')
        self.change_radii_all()
        if show_results:
            self.show_grid()
            self.show_grid('background')
            self.show_grid('cloud')
        self.minnaert_all(k=k)
        if show_results:
            self.map_all()
        else:
            self.map_all(showMaps = False)
    
    def write_img_to_png(self, filepath, wavlen, k = 0, cmap = 'bone', cmin=0, cmax=None):
        if cmax == None:
            cmax = 0
            for i in range(len(self.img_list)):
                tempmax = np.nanmax(self.img_list[i].img)
                if tempmax > cmax:
                    cmax = tempmax
        
        
        for elem in self.img_list:
            time = str(elem.expstart)
            time = time.replace('.','_')
            location = str(filepath)+time+'.png'
            plt.figure(figsize = (10,10))
            if k == 0:
                plt.imshow(elem.img, cmap=cmap)
                plt.text(elem.img.shape[1]*0.1,elem.img.shape[0]*0.1,str(elem.img_time.iso)+' - '+str(wavlen), fontsize=26, color='white')
            else:
                plt.imshow(elem.minnaert(k=k, img=elem.img), cmap=cmap)
                plt.text(elem.img.shape[1]*0.1,elem.img.shape[0]*0.1,str(elem.img_time.iso)+' - '+str(wavlen), fontsize=26)
            plt.clim(cmin,cmax)
            plt.axis('off')
            plt.savefig(location, bbox_inches='tight', pad_inches=0)
            plt.close()
            print(elem.identifier+' saved to png.')
    
    def write_to_fits(self, filepath, crop_of_img):
        for elem in self.img_list:
            x, y, z = img_to_xyz(elem.img.shape, elem.img.shape[0]/2, elem.img.shape[1]/2, elem.radius)
            theta, phi = xyz_to_thetaphi(x, y, z)
            theta_i, phi_i = xyz_to_thetaphi(x,y,z,theta=elem.sol_ang)
            latimg, longimg = xyz_to_longlat(x,y,z,180,elem.obslong, elem.obslat)
            mu = np.cos(theta*(np.pi/180))
            mu0 = np.cos(theta_i*(np.pi/180))
#             phase = np.zeros(mu.shape)
#             for i in range(mu.shape[0]):
#                 for j in range(mu.shape[1]):
#                     if mu[i,j] == np.nan:
#                         phase[i,j] = np.nan
#                     else:
#                         phase[i,j] = float(elem.eph['phi'])
            
            if crop_of_img == 'wide':
                hduimg = ast.PrimaryHDU(elem.img, header = elem.hdr)
                hdulong = ast.ImageHDU(longimg, header = elem.scihdr)
                hdulat = ast.ImageHDU(latimg)
                hdumu0 = ast.ImageHDU(mu0)
                hdumu = ast.ImageHDU(mu)
    #             hduphase = ast.ImageHDU(phase)
                hdulist = ast.HDUList([hduimg,hdulong,hdulat,hdumu,hdumu0])
                location = str(filepath)+'/'+elem.identifier+'_'+crop_of_img+'.fits'
                hdulist.writeto(location)
                print(elem.identifier+' images written to fits.')
            else:
                hduimg = ast.PrimaryHDU(elem.img, header = elem.hdr)
                hduback = ast.ImageHDU(elem.background, header = elem.scihdr)
                hducloud = ast.ImageHDU(elem.cloud_fill)
                hdulong = ast.ImageHDU(longimg)
                hdulat = ast.ImageHDU(latimg)
                hdumu0 = ast.ImageHDU(mu0)
                hdumu = ast.ImageHDU(mu)
    #             hduphase = ast.ImageHDU(phase)
                hdulist = ast.HDUList([hduimg,hduback,hducloud,hdulong,hdulat,hdumu,hdumu0])
                location = str(filepath)+'/'+elem.identifier+'_'+crop_of_img+'.fits'
                hdulist.writeto(location)
                print(elem.identifier+' images written to fits.')
            
        hdubackground = ast.PrimaryHDU(self.avg_background)
        hdulong = ast.ImageHDU(longimg)
        hdulat = ast.ImageHDU(latimg)
        hdumu0 = ast.ImageHDU(mu0)
        hdumu = ast.ImageHDU(mu)
        hdulist = ast.HDUList([hdubackground,hdulong,hdulat,hdumu,hdumu0])
        location = str(filepath)+'/avgbackground_'+crop_of_img+'.fits'
        hdulist.writeto(location)
        print('Average Background written to fits.')
            
    def write_map_to_fits(self, filepath):
        for elem in self.img_list:
            hduimg = ast.PrimaryHDU(np.flip(elem.mapped_img[0],0), header = elem.hdr)
            hduback = ast.ImageHDU(np.flip(elem.mapped_back[0],0))
            hducloud = ast.ImageHDU(np.flip(elem.mapped_cloud[0],0))
            hdulong = ast.ImageHDU(np.flip(elem.mapped_img[1],0))
            hdulat = ast.ImageHDU(np.flip(elem.mapped_img[2],0))
            hdulist = ast.HDUList([hduimg,hduback,hducloud,hdulong,hdulat])
            location = str(filepath)+elem.filter+'/'+elem.identifier+'_map.fits'
            hdulist.writeto(location)
            print(elem.identifier+' maps written to fits.')
            
    def avg_back_to_fits(self, filepath, filt):
        hduimg = ast.PrimaryHDU(self.avg_background[0])
        hdulist = ast.HDUList([hduimg])
        location = str(filepath)+str(filt)+'/average_background.fits'
        hdulist.writeto(location)
        print(elem.identifier+' maps written to fits.')
        
    def align_images(self, pad_value=0, centre_fac=10, **kwargs):
        args = []
        for elem in self.img_list:
            args.append(elem.img)
        img_list = []
        shape_list = []
        args = np.array(args)
        if len(args) == 1 and len(args[0].shape) > 2:
            args = args[0]

        for img in args:
            img = center_image(img, pad_value=pad_value, interp=centre_fac, **kwargs)
            img_list.append(img)
            shape_list.append(np.array(img.shape))

        shape_final = [max([s[0] for s in shape_list]),
                       max([s[1] for s in shape_list])]

        # Find needed offset for each image to align
        offset_list = []
        for shape in shape_list:
            offset_list.append(((shape_final - shape)/2).astype(int))

        output_list = []
        for idx in range(len(img_list)):
            img_aligned = np.zeros(shape_final, dtype=img_list[idx].dtype)
            img_aligned[:] = pad_value
            img_aligned[offset_list[idx][0]:offset_list[idx][0]+shape_list[idx][0],
                        offset_list[idx][1]:offset_list[idx][1]+shape_list[idx][1]] = img_list[idx]
            output_list.append(img_aligned)
        
        deinterp_list = []
        for img in output_list:
            deinterp_list.append(interp_image(img, 1/centre_fac))
        output_list = deinterp_list
            
        for i in range(len(output_list)):
            self.img_list[i].img = output_list[i]
            self.img_list[i].planicent = [output_list[i].shape[0]/2,output_list[i].shape[1]/2,self.img_list[i].radius]
     
    def despike_all(self, sigma = 1, times = 10):
        for elem in self.img_list:
            image = elem.img
            for i in range(times):
                image = exp_despike(image, sigma)
            elem.img = image
            print('image '+ elem.identifier +' has had hot pixels removed.')
    
    def set_planicentres_to_centre(self):
        for elem in self.img_list:
            elem.planicentre = [elem.img.shape[0]/2, elem.img.shape[1]/2, elem.radius]
            
print('image_list class initialised.')
