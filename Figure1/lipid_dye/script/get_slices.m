%     Copyright (C) 2024  Chia-Chou Wu
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.

%% initialize environment
clc; clear; close all

%% load images to slice
fname = "..\img\s41_cy5_ff1_10.tif";
% fname = "..\img\s41_ros_bgsub50.tif";
info = imfinfo(fname);
img = arrayfun(@(i)imread(fname,i,"Info",info),1:length(info),...
    'UniformOutput',0);

%% load setting for slicing images
oname = "..\data\setting_s41.mat";
load(oname)
mask = initiation.mask;
rot_angle = 225; % rotating for speed measurement

maskRotated = imrotate(mask,rot_angle); % rotated mask of initiation
S = regionprops(maskRotated,'Centroid');
center = fix(S.Centroid);

%% slice images
target_length = 4000; % um
px_length = ceil(target_length/1.266);

for j = 1:length(img)   

    imgRotated = imrotate(img{j},rot_angle);

    crop = imcrop(imgRotated,[center(1)-100,1,200,center(2)-1]);
    row_start = size(crop,1)-px_length+1;
    crop = crop(row_start:end,:,:);

    % save sliced images
    % imwrite(crop,"..\img\s41_cy5_ff1_10_sliced.tif","WriteMode","append")
    % imwrite(crop,"..\img\s41_ros_bgsub50_sliced.tif","WriteMode","append")
end
