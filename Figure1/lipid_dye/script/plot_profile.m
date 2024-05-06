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

%% load images to quantify
fname = "..\img\s41_ros_bgsub50.tif"; % lipid
info = imfinfo(fname);
img = arrayfun(@(i)imread(fname,i,"Info",info),1:length(info),...
    'UniformOutput',0);
dimg = arrayfun(@(i)img{i}-img{i-1},2:length(img),'UniformOutput',0);

fname = "..\img\s41_cy5_ff1_10.tif"; % cy5
info = imfinfo(fname);
img = arrayfun(@(i)imread(fname,i,"Info",info),1:length(info),...
    'UniformOutput',0);

oname = "..\data\setting_s41.mat";
load(oname)
mask = initiation.mask;

%% crop images
rot_angle = 270;
maskRotated = imrotate(mask,rot_angle);
S = regionprops(maskRotated,'Centroid');
center = fix(S.Centroid);

for j = 1:length(img)
    imgRotated = imrotate(img{j},rot_angle);
    crop = imcrop(imgRotated,[center(1)-100,16,200,2985]);
    imagelist_cy5{j} = crop;

    if j<length(img)
        imgRotated = imrotate(dimg{j},rot_angle);
        crop = imcrop(imgRotated,[center(1)-100,16,200,2985]);
        imagelist_dlp{j} = crop;
    end
end

%% plot mean intensity
figure('position',[680   558   560   200])
hold on
px_length = max(size(crop));
data = [];
color_sig = [0 126 255;254 212 57]./255;

for wSize = 40 % window size for moving average
    for j = 9 % frame to show
        colororder(color_sig)
        
        % cy5
        yyaxis("left")
        y = flipud(smoothdata(mean(imagelist_cy5{j},2),'movmean',wSize)) -...
            flipud(smoothdata(mean(imagelist_cy5{1},2),'movmean',wSize));
        y(y<0) = 0;
        y = normalize(y,1,'range');
        plot((1:px_length).*1.266,y,'color',color_sig(1,:),'linewidth',2)
        ylim([-0.05 1.05])        
        data = [(1:px_length)'.*1.266,y];
        
        % ros
        yyaxis('right')
        y = flipud(smoothdata(mean(imagelist_dlp{j-1},2),'movmean',wSize));
        y(y<0) = 0;
        y = normalize(y,1,'range');
        plot((1:px_length).*1.266,y,'color',color_sig(2,:),'linewidth',2)
        ylim([-0.05 1.05])
        data = [data,y];

        xlim([0 length(y)*1.266])
        xticks(0:1000:3000)
        xticklabels(0:3)
        xlabel("Distance (mm)")
 
    end
    hold off    
end
set(gca,'FontSize',16,'LineWidth',1,'TickDir','out','YColor','k')

