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

fname = ["DFO_24_s7" "DFO_24_s14" "FAC_24_s1" "FAC_24_s5",...          % 1 2 3 4
    "LY294002_24_s19" "LY294002_24_s21",...                            % 5 6
    "Dasatinib_20211113_s24" "Dasatinib_20211117_s4",...               % 7 8
    "GKT_24_20210827_s9" "GKT_24_20210827_s14" "GKT_24_20211113_s7",...% 9 10 11
    "GKT_24_20211113_s8" "GKT_24_20211117_s11"]+"_c1_ff1_10";               % 12 13
th = [1 18e-3 1 10e-3,...
      1 1,...
      1 1,...
      1 1     1,...
      1 1];
frame1_treat = [12 12 16 16,...
                13 13,...
                13 13,...
                13 13 13,...
                13 13];
for well = 2
    % cy5
    imgname = "..\img\"+fname(well)+".tif"; sig = "cy5"; % ros
    info = imfinfo(imgname);
    img = arrayfun(@(i)imread(imgname, i, 'Info', info),1:length(info),...
        'UniformOutput',0);
    
    % create make of death cell population
    mask = cellfun(@(x)imbinarize(x,th(well)),img,'UniformOutput',0);
    mask = cellfun(@(x)bwareafilt(x,[0 1000]),mask,'UniformOutput',0);
    mask = cellfun(@(x)imdilate(x,strel('disk',25)),mask,'UniformOutput',0);
    mask = cellfun(@(x)bwareaopen(x,50000),mask,'UniformOutput',0);
    mask = cellfun(@(x)~bwareaopen(~x,5000),mask,'UniformOutput',0);

    % smooth outlines
    diskRadius = 60;
    se = strel('disk', diskRadius, 0);
    kernel = se.Neighborhood / sum(se.Neighborhood(:));
    blurredImage = cellfun(@(i)conv2(double(i), kernel, 'same'),mask,'UniformOutput',0);
    bo = cellfun(@(i) bwboundaries(i > 0.3), blurredImage,'UniformOutput',0);
    bdy = cell(1,length(bo));
    for i = 1:length(bo)
        if ~isempty(bo{i})
            bdy{i} = bo{i}{1};
            if well==4 && i==17
                bdy{i} = bo{i}{2};
            end
        end
    end
    toc
    
    % identify the photoinduction site
    i = find(cellfun(@(i)~isempty(i),bdy),1);
    BW = poly2mask(bdy{i}(:,2),bdy{i}(:,1),size(img{1},1),size(img{1},2));
    centroid = regionprops(BW,'Centroid');
    centroid = fix(centroid.Centroid);    
    coord_intersect = [];
    
    % plot outlines
    close
    figure
    imshow(img{frame1_treat(well)},[10 1500],'border','tight');
    for i = 1:size(bdy,2)
        hold on
        if ~isempty(bdy{i})
            x = bdy{i}(:,2);
            y = bdy{i}(:,1);
            isSetToNan = (x==1|x==size(img{1},2))|(y==1|y==size(img{1},1));
            x(isSetToNan) = nan;
            y(isSetToNan) = nan;
            if i == frame1_treat(well)
                plot(x,y,'m','LineWidth',2)
            else
                plot(x,y,'w','LineWidth',2)
            end
            [xi,yi,ii] = polyxpoly([1 5000],[4621 4621],x,y);
            if ~isempty(ii)
                [~,ind] = max(xi);
                coord_intersect = [coord_intersect;i xi(ind) ii(ind,2)];
            end
        end
        hold off
    end
    
    viscircles(centroid,197,'Color','c','EnhanceVisibility',0);
    set(gca,'position',[0 0 1 1],'units','normalized')
    axis ij

    hold on
    plot(coord_intersect(:,2),4621.*ones(size(coord_intersect,1)),'ys-')
    hold off
    
end