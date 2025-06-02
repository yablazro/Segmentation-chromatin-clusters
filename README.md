# Segmentation-chromatin-clusters
Segmentation algorithm to identify chromatin clusters in fluorescent microscope images.

READING INTENSITY IMAGES FOR CHROMATIN SEGMENTATION           
---------------------------------------------------
Author: María del Valle Blázquez-Romero
Contact info: yablazro@gmail.com
Date: 17 Feb. 2025
Research group: M2BE - I3A - Universidad de Zaragoza

SUMMARY

This code reads image files from a folder, converts them into intensity
images (double) ranging from 0 to 1, and process them through the user
function called "SegmentationRegionsSimulations.m", which segments bright
regions within the image by applying a global threshold together with a 
locally adaptive threshold.

STEPS
1. Prepare a folder containing all the individual images that you want to
   process.
2. Press run and introduce the image file type in the input window that 
   will appear (make sure that all your files are of the same type).
3. In the next window, enter the value for the parameters defining the
   global and local segmentation of bright regions and the Gaussian Blur
   filter: maxT, minT, sens and devstd (see definitions below).
4. Select the normalization method (explained below) to post-process the 
   images. Generally use individual pixel normalization.
5. Indicate the input directory by selecting the folder containing the
   images to analyze.
6. The results will appear in different figures:
   - Figure 1, raw images.
   - Figure 2, normalized images.
   - Figure 3, Binary mask for bright regions segmentation.
   - Figure 4, Binary mask for total nucleus area segmentation (assumed
     to be a circunscribed circle).
   - Figure 5, Final segmentation result, overlapping the contours of
     bright regions on the normalized image.

SEGMENTATION PARAMETERS

Segmentation parameters (maxT, minT and sens) may need to be optimized
for different images. Gaussian Blur filter can be adjusted with variable 
devstd.

maxT   => Threshold to include highest-intensity pixels in global 
          segmentation.
          
minT   => Threshold to ignore lowest-intensity pixels in local
          segmentation even if the filters include them.
          
sens   => Sensitivity parameter defined for function "adaptthresh" to
          segment high-intensity regions locally.
          
devstd => Specifies the standard deviation used for the gaussian blur
          fliter ("imgaussfilt").

NORMALIZATION METHODS

Three option to define the method used to normalize the raw image:

1 => Individual pixel normalization (I_norm=(I-Imin)/Imean), where I is 
     the raw intensity of the pixel, Imin minimum intensity within the 
     frame and Imean the mean intensity of the frame ignoring zeros.
     
2 => MATLAB function "imadjust"

3 => None

- Go to "SegmentationRegionsSimulations.m" for the segmentation code.
