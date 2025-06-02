%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   CODE FOR HIHGLY-CONDENSED CHROMATIN SEGMENTATION IN INTENSITY IMAGES  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Author: María del Valle Blázquez-Romero
% Contact info: yablazro@gmail.com
% Date: 17 Feb. 2025
% Research group: M2BE - I3A - Universidad de Zaragoza

% -------------------------------------------------------------------------
% INPUT VARIABLES
% -------------------------------------------------------------------------
% raw        => matrix array containing the raw intensities for the pixels
%               within the nucleus for ONE image.
% n_imgs     => number of individual images to segment to keep track.
% norm_indx  => defines the method used to normalize the raw image:
%               1 = Individual pixel normalization (I_norm=(I-Imin)/Imean)
%               2 = MATLAB function "imadjust"
%               3 = None
% maxT       => Threshold to include highest-intensity pixels in global 
%               segmentation.
% minT       => Threshold to ignore lowest-intensity pixels in local
%               segmentation even if the filters include them.
% sens       => Sensitivity parameter defined for function "adaptthresh" to
%               segment high-intensity regions locally.
% devstd     => Specifies the standard deviation used for the gaussian blur
%               fliter ("imgaussfilt").

% -------------------------------------------------------------------------
% OUTPUT VARIABLES
% -------------------------------------------------------------------------
% imgNorm_final => matrix array containing the normalized intensities for 
%                  the pixels within the nucleus for ONE image.
% segNucleus    => mask as a binary matrix array with the result for the
%                  total nucleus area segmentation.
% segHet        => mask as a binary matrix array with the result for the
%                  segmentation of highly-condensed chromatin regions.
% area          => structure containing the numerical result in #pixels for
%                  total nucleus area (A_tot) and total highly-condensed
%                  area (A_con)
% intNorm       => structure containing the numerical result as "double"
%                  for mean intensity and maximum intensity after
%                  normalization.
% morph         => structure containing morphological parameters for each
%                  of the segmented highly-condensed chromatin regions.

% -------------------------------------------------------------------------
% Segmentation parameters (maxT, minT and sens) may need to be optimized
% for different images.
% -------------------------------------------------------------------------

function [imgNorm_final,segNucleus,segHet,area,intNorm,morph] = SegmentationRegionsSimulations(img,n_imgs,norm_indx,maxT,minT,sens,devstd)

% Measure the intensities within the frame
intMean = mean(img(img~=0));                                                % Does NOT consider zeros in the background
intMin = min(img(img~=0),[],'all','omitnan');                               % Does NOT consider zeros in the background

fprintf('Mean intensity w/o 0 for image #%d = %4.4f \n',n_imgs,intMean)
fprintf('Minimum intensity w/o 0 for image #%d = %4.4f \n\n',n_imgs,intMin);

switch norm_indx
    case 1
        % Individual pixel normalization: I_norm = (I-Imin)/Imean
        % where I is the raw intensity of the pixel, Imin minimum intensity within
        % the frame and Imean the mean intensity of the frame
        imgNorm = (img - intMin)./intMean;
    case 2
        % Adjust image contrast
        imgNorm = imadjust(img);
    case 3
        imgNorm = img;
end

% Gauss filter to smooth edges
imgNorm_gauss = imgaussfilt(imgNorm,devstd);

% Segmentation of the nucleus area by drawing a circunscribed circle
dim = size(img);
area.Nucleus = pi*(dim(1)/2)^2;
r = dim(1)/2;                                                               % Radius of the circunscribed circle
x = r;                                                                      % X coordinate for the center of the circle
y = r;                                                                      % Y coordinate for the center of the circle
th = 0:pi/50:2*pi;
x_circle = r * cos(th) + x;
y_circle = r * sin(th) + y;
segNucleus = poly2mask(x_circle,y_circle,dim(1),dim(2));                    % Create binary mask for the circunscribed circle defining the nucleus

% Segmentation of highly-condensed chromatin regions
threshold_hetlocal = adaptthresh(imgNorm_gauss,sens);                       % Define adaptative threshold to segment locally high-intensity regions
threshold_hetlocal(imgNorm_gauss <= minT) = 1;                              % Sets the threshold to 1 to ignore low-intensity pixels in the adaptative filter

imgNorm_final = imgNorm_gauss.*segNucleus;                                  % Outside of the nucleus = 0

segHetLocal = imbinarize(imgNorm_final, threshold_hetlocal);                % Create mask for local segmentation
segHetGlobal = imbinarize(imgNorm_final, maxT);                             % Create mask for global segmentation (considers only pixels greater than max_thres)

segHet = segHetLocal + segHetGlobal;                                        % Sum local and global masks to obtain the final segmentation


% Calculate areas
area.Nucleus = nnz(segNucleus);                                             % # pixels different from zero in the nucleus converted to um^2
area.Het = nnz(segHet);                                                     % # pixels different from zero after region segmentation converted to um^2 

% Measure normalized intensities
intNorm.Mean = mean(imgNorm_final);
intNorm.Max = max(imgNorm_final,[],'all');

% Measure properties of connected regions (individual highly-condensed
% chromatin regions)
CC = bwconncomp(segHet);
morph = regionprops(CC,"all");
end