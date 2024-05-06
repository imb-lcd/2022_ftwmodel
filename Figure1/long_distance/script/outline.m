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

% load cyanHot colormap
load("cyanHot.mat")

% load image
imgname = "../img/s19_cy5_ff1_10.tif";
info = imfinfo(imgname);
imgs = arrayfun(@(s)imread(imgname,s,"Info",info),3:19,'UniformOutput',0);

% mask death population
th = [535 535 535 535 555 555 555 555 575 575 575 575,570 500 500 495 490];
for i = 1:length(imgs)

    th_high = th(i);
    th_area_lb = 0;
    th_area_ub = 2000;
    r_dilation = 25;
    r_erosion = 5;

    % binarize the image
    mask = imbinarize(imgs{i},th_high/2^16);
    mask = bwareafilt(mask,[th_area_lb th_area_ub]);

    % dilation
    maskD = bwdist(mask)<r_dilation;
    maskD = ~bwareaopen(~maskD,th_area_ub);

    % erosion
    maskDS = bwdist(~maskD)>r_erosion;

    % size filtering
    maskDS = bwareaopen(maskDS,50000);
    maskDS = ~bwareaopen(~maskDS,50000);    

    % smooth outlines
    diskRadius = min([20+(i-1)*10 60]);
    se = strel('disk', diskRadius, 0);
    kernel = se.Neighborhood / sum(se.Neighborhood(:));
    blurredImage = cellfun(@(i)conv2(double(i), kernel, 'same'),{maskDS},'UniformOutput',0);
    maskDS_smoothed = cellfun(@(x)imfill(x > 0.2,'holes'),blurredImage,'UniformOutput',0);
    maskDS = maskDS_smoothed{1};
    
    % keep the biggest region
    L = bwlabel(maskDS);
    value = setdiff(min(L(:)):max(L(:)),0);
    area = arrayfun(@(v)nnz(L==v),value);
    [~,index] = max(area);
    maskDS = L==value(index);

    MASK{i} = maskDS;
end

% create landscape
landscape = uint8(cat(3,MASK{:}));
[a,Z]=max(landscape,[],3);
Z(a==0) = size(landscape,3)+1;

% plot outlines on the image
close all
figure
imshow(medfilt2(imgs{10},[2 2]),[235 235+255],'Border','tight')
colormap(cyanHot)
hold on
for i = 1:max(Z(:))
    mask_ros = Z<=i;
    if any(mask_ros(:))
        outline_ros = bwboundaries(mask_ros);
        for j = 1:length(outline_ros)
            bdy = outline_ros{j};
            bdy(sum(ismember(bdy,[1 5000]),2)>0,:) = nan;
            if i==10
                plot(bdy(:,2),bdy(:,1),'y','linewidth',2)
            else
                plot(bdy(:,2),bdy(:,1),'w','linewidth',2)
            end
        end
        
    end
end
xlim([1 4500])
ylim([501 5000])

S = regionprops(Z<=1,'centroid');
center = fix(S.Centroid);
viscircles(center,197,'Color','r','EnhanceVisibility',false)
pgon = polyshape(center(1)-450+[0 900 900 0],...
    center(2)-450+[0 0 900 900]);
hold on
p = plot(pgon,'FaceColor','none','EdgeColor','y','LineWidth',2,'LineStyle','--');
hold off
% exportgraphics(gcf,"outline_global.pdf","ContentType","vector");

xlim(center(1)-450+[0 900])
ylim(center(2)-450+[0 900])
% exportgraphics(gcf,"outline_local.pdf","ContentType","vector");