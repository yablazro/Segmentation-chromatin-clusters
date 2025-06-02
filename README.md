# Segmentation-chromatin-clusters
Segmentation algorithm to identify chromatin clusters in fluorescent microscope images. The main code to run the segmentation is called "ReadImg.m" while "SegmentationChromatinClusters.m" corresponds to a user function containing the segmentation algorithm. Find below a summary for the work flow of each code. The details and step-by-step process for each code is contained as a comment within the files ".m".

Author: María del Valle Blázquez-Romero

Contact info: yablazro@gmail.com

Date: 17 Feb. 2025

Research group: M2BE - I3A - Universidad de Zaragoza


ReadImg.m           
----------
This code reads image files from a folder, converts them into intensity
images (double) ranging from 0 to 1, and process them through the user
function called "SegmentationRegionsSimulations.m".

SegmentationChromatinClusters.m           
--------------------------------
This code segments bright regions within the fluorescent images by applying a global threshold together with a 
locally adaptive threshold.
