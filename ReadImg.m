%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           READING INTENSITY IMAGES FOR CHROMATIN SEGMENTATION           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Author: María del Valle Blázquez-Romero
% Contact info: yablazro@gmail.com
% Date: 17 Feb. 2025
% Research group: M2BE - I3A - Universidad de Zaragoza

% -------------------------------------------------------------------------
% SUMMARY
% -------------------------------------------------------------------------
% This code reads image files from a folder, converts them into intensity
% images (double) ranging from 0 to 1, and process them through the user
% function called "SegmentationRegionsSimulations.m", which segments bright
% regions within the image by applying a global threshold together with a 
% locally adaptive threshold.

% -------------------------------------------------------------------------
% STEPS
% -------------------------------------------------------------------------
% 1. Prepare a folder containing all the individual images that you want to
%    process.
% 2. Press run and introduce the image file type in the input window that 
%    will appear (make sure that all your files are of the same type).
% 3. In the next window, enter the value for the parameters defining the
%    global and local segmentation of bright regions and the Gaussian Blur
%    filter: maxT, minT, sens and devstd (see definitions below).
% 4. Select the normalization method (explained below) to post-process the 
%    images. Generally use individual pixel normalization.
% 5. Indicate the input directory by selecting the folder containing the
%    images to analyze.
% 6. The results will appear in different figures:
%    - Figure 1, raw images.
%    - Figure 2, normalized images.
%    - Figure 3, Binary mask for bright regions segmentation.
%    - Figure 4, Binary mask for total nucleus area segmentation (assumed
%      to be a circunscribed circle).
%    - Figure 5, Final segmentation result, overlapping the contours of
%      bright regions on the normalized image.

% -------------------------------------------------------------------------
% SEGMENTATION PARAMETERS
% -------------------------------------------------------------------------
% Segmentation parameters (maxT, minT and sens) may need to be optimized
% for different images. Gaussian Blur filter can be adjusted with variable 
% devstd.
%
% maxT   => Threshold to include highest-intensity pixels in global 
%           segmentation.
% minT   => Threshold to ignore lowest-intensity pixels in local
%           segmentation even if the filters include them.
% sens   => Sensitivity parameter defined for function "adaptthresh" to
%           segment high-intensity regions locally.
% devstd => Specifies the standard deviation used for the gaussian blur
%           fliter ("imgaussfilt").

% -------------------------------------------------------------------------
% NORMALIZATION METHODS
% -------------------------------------------------------------------------
% Three option to define the method used to normalize the raw image:
% 1 => Individual pixel normalization (I_norm=(I-Imin)/Imean), where I is 
%      the raw intensity of the pixel, Imin minimum intensity within the 
%      frame and Imean the mean intensity of the frame ignoring zeros.
% 2 => MATLAB function "imadjust"
% 3 => None

% -------------------------------------------------------------------------
% Go to "SegmentationRegionsSimulations.m" for the segmentation code.
% -------------------------------------------------------------------------

close all;
clear;
clc;

prompt = {'Enter file type (png,jpg,tif...):'};
dlgtitle = 'Data configuration';
fieldsize = [1 45];
definput = {'png'};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);

filetype = answer{1,1};

prompt = {'Enter maximum threshold to consider high-intensity pixels:',...
    'Enter minimum threshold to ignore low-intensity pixels:',...
    'Enter sensitivity parameter for local segmentation with "adptthresh":'...
    'Enter standard deviation for Gaussian Blur Filter:'};
dlgtitle = 'Segmentation parameters';
fieldsize = [1 45; 1 45; 1 45; 1 45];
definput = {'0.90','0.60','0.50','2.5'};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);

maxT = str2num(answer{1,1});
minT = str2num(answer{2,1});
sens = str2num(answer{3,1});
devstd = str2num(answer{4,1});

norm_method = {'Indv. pixel norm','Imadjust','None'};
[norm_indx,~] = listdlg('PromptString',{'Select the normalization method to be used',''},...
    'SelectionMode','single','ListString',norm_method);

folder = uigetdir('Select the path containing the results from Matlab analysis');
filetype = ['*.',filetype];
files = dir(fullfile(folder,filetype));
units = 'pix';

imgs_raw = cell(length(files),1);
imgs_norm = cell(length(files),1);
segNuclei = cell(length(files),1);
segHets = cell(length(files),1);
areasNuclei = cell(length(files),1);
areasHet = cell(length(files),1);
areasRel = cell(length(files),1);
intsMean = cell(length(files),1);
intsMax = cell(length(files),1);

figure(1);
sgtitle('Raw Images')
figure(2);
sgtitle('Normalized Images')
figure(3);
sgtitle('Bright regions segmentation')
figure(4)
sgtitle('Nucleus area segmentation');
figure(5)
sgtitle(['Segmentation result with maxT = ',num2str(maxT),', minT = ',num2str(minT),', sens = ',num2str(sens),' and devstd = ',num2str(devstd)]);


for f=1:length(files)
    raw = imread([folder '/' files(f).name]);
    imgs_raw{f,1} = im2double(raw(:,:,1));
    
    [imgNorm_final,segNucleus,segHet,area,intNorm,morph] = SegmentationRegionsSimulations(imgs_raw{f,1},f,norm_indx,maxT,minT,sens,devstd);
    
    imgs_norm{f,1} = imgNorm_final;
    segNuclei{f,1} = segNucleus;
    segHets{f,1} = segHet;
    areasNuclei{f,1} = area.Nucleus;
    areasHet{f,1} = area.Het;
    intsMean{f,1} = intNorm.Mean;
    intsMax{f,1} = intNorm.Max;
    
    n_plot = ceil(sqrt(length(files)));
    m_plot = n_plot;
    
    areasRel{f,1} = areasHet{f,1}/areasNuclei{f,1}*100;
    
    figure(1)
    hold on;
    subplot(n_plot,m_plot,f)
    imshow(imgs_raw{f,1})
    hold on;
    
    figure(2)
    hold on;
    subplot(n_plot,m_plot,f)
    imshow(imgs_norm{f,1})
    hold on;
    
    figure(3)
    hold on;
    subplot(n_plot,m_plot,f)
    imshow(segHets{f,1})
    hold on;
    title(['A_{het} = ',num2str(areasHet{f,1}),' ',units]);

    figure(4)
    hold on;
    subplot(n_plot,m_plot,f)
    imshow(segNuclei{f,1})
    hold on;
    title(['A_{tot} = ',num2str(areasNuclei{f,1}),' ',units]);
    
    figure(5)
    hold on;
    subplot(n_plot,m_plot,f)
    imshow(imgs_norm{f,1})
    hold on;
    visboundaries(segNuclei{f,1},'Color',[1 0 0])
    hold on;
    visboundaries(segHets{f,1},'Color',[0 1 0])
    hold on;
    title(['A_{het}/A_{tot} = ',num2str(areasRel{f,1}),'%']);
    
    clear raw nuclei_seg gImage nuclei_contours customColormap mapHet rgbHet color NulceiHet
end

fprintf('PARAMETERS FOR BRIGHT REGION SEGMENTATION:')
fprintf('\nMaximum threshold (maxT) = %2.2f',maxT)
fprintf('\nMinimum threshold (minT) = %2.2f',minT)
fprintf('\nSensitivity (sens) = %2.2f',sens)
fprintf('\n\nGAUSSIAN BLUR:')
fprintf('\nStandard deviation (devstd) = %2.2f\n\n',devstd)