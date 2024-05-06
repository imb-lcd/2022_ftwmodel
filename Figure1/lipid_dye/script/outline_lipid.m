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

load("cyanHot.mat")

isSave = 0;
well = "41";
list = dir("..\img\s"+well+"\for_outline\ros*.tif");

%%
MASK = cell(size(list));

for i = 1:length(list)
    close all

    % load image for outlining wave
    fname = fullfile(list(i).folder,list(i).name);
    img = flip(flip(imread(fname),2),1);

    % if the outline exists, load the setting; otherwise use the default
    if exist("..\data\mask_s"+well+".mat","file")
        load("..\data\mask_s"+well+".mat","th_high","th_area_lb","th_area_ub","r_dilation","r_erosion")
        isSave = 0;
    else
        % ros
        th_high = 20000; % 18000 592
        th_area_lb = 20;
        th_area_ub = 1000;
        r_dilation = 35;
        r_erosion = 25;
    end

    % binarize the image
    mask = imbinarize(img,th_high/2^16);
    mask = bwareafilt(mask,[th_area_lb th_area_ub]);
    figure("Name","High threshold"); imshow(mask,'Border','tight')

    % dilation
    %     maskD = imfill(bwdist(mask)<r_dilation,'holes');
    maskD = bwdist(mask)<r_dilation;
    maskD = ~bwareaopen(~maskD,th_area_ub);

    % erosion
    maskDS = bwdist(~maskD)>r_erosion;

    % size filtering
    maskDS = bwareaopen(maskDS,5000);
    maskDS = ~bwareaopen(~maskDS,100000);

    % smooth outlines
    H = fspecial('average',60);
    blurredImage = imfilter((maskDS),H,'replicate');
    maskDS_smoothed = imfill(blurredImage > 0.2,'holes');
    maskDS = maskDS_smoothed;

    % smooth outlines
    diskRadius = min([20+(i-1)*10 60]); % less smooth
    se = strel('disk', diskRadius, 0);
    kernel = se.Neighborhood / sum(se.Neighborhood(:));
    blurredImage = cellfun(@(i)conv2(double(i), kernel, 'same'),{maskDS},'UniformOutput',0);
    maskDS_smoothed = cellfun(@(x)imfill(x > 0.2,'holes'),blurredImage,'UniformOutput',0);
    maskDS = maskDS_smoothed{1};

    % keep the biggest region
    L = bwlabel(maskDS);
    value = setdiff(min(L(:)):max(L(:)),0);
    if isempty(value)
        continue
    end
    area = arrayfun(@(v)nnz(L==v),value);
    [~,index] = max(area);
    maskDS = L==value(index);

    MASK{i} = maskDS;
end

%% show landscape
landscape = uint8(cat(3,MASK{:}));
[a,Z]=max(landscape,[],3);
Z(a==0) = size(landscape,3)+1;

close all
imshow(Z,[min(Z(:)) max(Z(:))],'border','tight')

%%
i = 7
fname = fullfile(list(i).folder,list(i).name);
img = flip(flip(imread(fname),2),1);

close all
figure
imshow(imadjust(img),'Border','tight')
% imshow(img,[212 467],'Border','tight')
colormap(cyanHot)
hold on
for i = 1:max(Z(:))
    mask_ros = Z<=i;
    if any(mask_ros(:))
        outline_ros = bwboundaries(mask_ros);
        for j = 1:length(outline_ros)
            bdy = outline_ros{j};
            bdy(sum(ismember(bdy,[1 3000]),2)>0,:) = nan;

            plot(bdy(:,2),bdy(:,1),'y','linewidth',2)
        end
        
    end
end

S = regionprops(Z<=1,'centroid');
center = fix(S.Centroid);

viscircles(center,197,'Color','r','EnhanceVisibility',0,'LineStyle','--')

pgon = polyshape([647 647+1499 647+1499 647],...
    [971 971 971+1499 971+1499]);
hold on
p = plot(pgon,'FaceColor','none','EdgeColor',[.7 .7 .7],'LineWidth',2,'LineStyle','--');
hold off
% print('-painters',gcf,"..\fig\outline_s41_ros_global.svg",'-dsvg')

xlim([647 647+1499])
ylim([971 971+1499])
% print('-painters',gcf,"..\fig\outline_s41_cy5_local.svg",'-dsvg')