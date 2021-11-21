---Pipeline Construction Guide---

The code in the image classes here was designed to allow for easy construction of a processing pipeline, which turns a HST or HLA image of Neptune from the MAST archive (https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html) into an "archival" data format. 
This file shows the intended pipeline design, for the user to form a pipeline for use on their own system, which is likely to have a different file structure. First the "recommended" file structure is given, then the two versions of the image processing pipeline which I have created. One for a minimum crop image, and one for a "close" crop of the image, to an "image radius" of 2 Neptune Radii.

---Recommended File Structure---

Planet_Image_Processing (directory to contain all code and data)
|- ancilliary
|           \- filter_fluxes (a directory created by Dr Michael Roman, containing the HST WFPC2 and WFC3 filter passthrough fluxes for I/F calculation purposes)
|                          |- filter_solar_fluxes.txt
|                          \- README.txt
|           
|
|- data (Store the HST and HLA pipeline data here, for the code to access, in whatever format you choose. My preferred style of file structure is shown below)
|     \- IMAGE WAVELENGTHS, one folder per FILTER CENTRAL WAVELENGTH, in format ***nm
|             \- *MAST IMAGE IDENTIFIER - ALL FOLDERS FROM MAST* (many folders here, each one containing a single image)
|                    \- *MAST IMAGE IDENTIFIER*.fits
|
|- processed_data
|               |- fits
|               |     \- IMAGE WAVELENGTHS, one folder per FILTER CENTRAL WAVELENGTH, in format ***nm
|               |                                                                                   \- IDENTIFIER_CROP.fits
|               |- gifs
|               |     \- (store any gifs created from the images here, gifs generated from png images)
|               |
|               \- pngs (filenames are the modified julian date of the image, to allow for easy and automatic file ordering)
|                     |- imgs
|                     |     \- IMAGE WAVELENGTHS, one folder per FILTER CENTRAL WAVELENGTH, in format ***nm
|                     |                                                                                   \- *****_*********.png
|                     |- background
|                     |           \- IMAGE WAVELENGTHS, one folder per FILTER CENTRAL WAVELENGTH, in format ***nm
|                     |                                                                                   \- *****_*********.png
|                     \- cloud_fill 
|                                 \- IMAGE WAVELENGTHS, one folder per FILTER CENTRAL WAVELENGTH, in format ***nm
|                                                                                                         \- *****_*********.png
|
\- processing_code
                 |- func
                 |     |- image_class.py
                 |     |- image_list_class.py
                 |     \- MapFunc.py
                 |
                 |- tools
                 |      \- OLIVER KING's ASTROTOOLS MODULE, ACCESS AT https://github.com/ortk95/astro-tools
                 |
                 \- *******.ipynb OR *******.py (pipeline files of your own creation)
                 
                 
---Recommended Pipeline Structure---
There are 2 recommended pipeline structures, one for "wide" (minimum loss) cropping and one for "close" (2*Neptune Radius) cropping. This is due to inaccuracy with image centring on the wide cropped image, which means that the background and cloud separation routine doesn't give accurate results. Therefore, I've removed that section of the pipeline for the wide crop image, choosing to replace the "background" and "cloud_fill" image extensions with the original image in that case.

--- Close Crop Pipeline --- python code

from func.image_list_class import *
files = filesearch('(YOUR FILEPATH HERE)/Planet_Image_Processing/data/(WAVELENGTH)nm/*/*_drz.fits')
bad = [] #after running through the pipeline once, and determining the images which aren't applicable to use for data, use this list to list the identifiers of the bad images as strings
img = image_list()
img.add_images(files)
img.trim_edge_all()
img.img_goodness(badimgs = bad)
img.despike_all()
img.calib_all()
img.planicent_all()
img.shift_all()
img.normalise_radius()
img.crop_all(crop='close')
img.rotate_all()
img.align_images()
img.set_planicentres_to_centre()
img.crop_all(crop='close')
img.isolate_backgrounds(cloudAlbedo = 0.9)
img.write_to_fits(filepath = 'C:/Users/jamie/Documents/Planet_Image_Processing/processed_data/fits/(WAVELENGTHnm)/', crop_of_img='close')

--- Wide Crop Pipeline --- python code

from func.image_list_class import *
files = filesearch('(YOUR FILEPATH HERE)/Planet_Image_Processing/data/(WAVELENGTH)nm/*/*_drz.fits')
bad = [] #after running through the pipeline once, and determining the images which aren't applicable to use for data, use this list to list the identifiers of the bad images as strings
img = image_list()
img.add_images(files)
img.trim_edge_all()
img.img_goodness(badimgs = bad)
img.despike_all()
img.calib_all()
img.planicent_all()
img.shift_all()
img.normalise_radius()
img.crop_all(crop='wide')
img.rotate_all()
img.write_to_fits(filepath = 'C:/Users/jamie/Documents/Planet_Image_Processing/processed_data/fits/(WAVELENGTHnm)/', crop_of_img='wide')
