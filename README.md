# Planet_Image_Processing
This repository contains code for processing and mapping images of planetary bodies. Designed for images of neptune taken by Hubble Wide Field Instruments, however could likely be adapted to fit other telescopes.

Dependancies:
- Oliver King's astro-tools (https://github.com/ortk95/astro-tools) and all dependancies
- numpy
- scipy
- matplotlib
- astropy
- photutils

In it's current state, this code takes an image, and turns it into an archival fits format. Code to map an image is present, however it is not recommended to be used, I am planning on creating a "archive_image" class, which opens the archival format, and can perform photometry, and map a given image. I am also planning on creating an "archive_image_list" class, to perform these steps for a series of images, and also make note of groups of images within the dataset with close image taken timings, for creation of a "Neptune Rotation", allowing for full longitude mapping range.

Before using any of the code, it is recommended to also read pipeline_construction_readme.txt
