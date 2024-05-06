%% initialize environment
clc; clear; close all

% set the filename of the output
oname = "..\data\setting_s19.mat";

% default setting
% frame_index: frame of initiation
% signal: which signal is going to be used for identify initiation
% range_area: filter out the region whose area is outside the range
% r_dilation/r_erosion: connect dead cells
% th_mask: the area of connected regions smaller than th_mask will
%          be removed.
% boundary_index: the index of the region where the initiation
%                 locates
initiation = struct('frame_index',3,'signal',"cy5",'th_intensity',0.9,...
    'range_area',[20 1000],'r_dilation',15,'r_erosion',10, ...
    'th_mask',5000,'boundary_index',1);

% load the image with the initiation to be identify
imgname = "..\img\s19_cy5_ff1_10.tif";
img = imread(imgname, initiation.frame_index, 'Info', imfinfo(imgname));

% binarize the image
im_adjusted = imadjust(img,stretchlim(img));
mask = imbinarize(im_adjusted,initiation.th_intensity);

% filter out the regions with abnormal area
mask = bwareafilt(mask,initiation.range_area);
% imshow(BW)

% dilate the threholding result to connect dead cells
maskD = imfill(bwdist(mask)<initiation.r_dilation,'holes');
maskDS = bwdist(~maskD)>initiation.r_erosion;

% remove small objects after connection
maskDS = bwareaopen(maskDS,initiation.th_mask);

% overlay the remaining regions with the image
L = bwlabel(maskDS);
B = labeloverlay(imadjust(medfilt2(img,[3 3]),[234 489]./2^16),bwlabel(maskDS));

% check initiation
figure('WindowState','maximized')
tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
nexttile
imshow(L,[0 max(L(:))],'Border','tight');
title('Label matrix')
nexttile
imshow(B)
title('Image overlaid with the label matrix')

% mark the initiation center
maskDS = L==initiation.boundary_index;
hold on
stats = regionprops(maskDS,'Centroid');
plot(stats(1).Centroid(1),stats(1).Centroid(2),'rx','MarkerSize',20)

% save mask into setting
initiation.mask = maskDS;

save(oname,"initiation")
