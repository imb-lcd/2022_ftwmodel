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

tbl = readtable("..\data\parameters.csv");
line_alpha = [0.15 0.2 0.3 0.5 0.6 0.7 1]; % transparency
line_color = [repmat([255 209 0]./255,7,1),line_alpha']; % color

% load gap size
load("..\data\gapSize.mat")

%% for each well, get the intensity profile of ROS
data = [];
for i = [1 2 7 8]
    
    % load processed images
    fname = "..\img\s"+tbl.wells(i)+"_ros_ff1_10_aligned_cropped.tif"; % ros
    info = imfinfo(fname);
    ros = arrayfun(@(i)imread(fname,i,"Info",info),1:length(info), ...
        'UniformOutput',0);
    
    % deal with the black frames in the well 89 and 30
    if tbl.wells(i)==89
        ros([4 8]) = [];
    end
    if tbl.wells(i)==30
        ros([6 10]) = [];
    end

    % get differential ROS
    ros = [ros(1),ros];
    dros = arrayfun(@(i)ros{i}-ros{i-1},2:length(ros),'UniformOutput',0);
    
    % get the local view
    dros = dros(tbl.frame_strt_global(i):tbl.frame_stop_global(i));
    dros = dros(tbl.frame_strt_local(i):tbl.frame_stop_local(i));
    dros = cellfun(@(x) imcrop(imrotate(x,abs(tbl.angle(i)),'nearest','crop'),...
        [tbl.x_local(i),tbl.y_local(i),200,round(2000/tbl.pxSize(i))-1]),dros,'UniformOutput',0);
    
    % get intensity profile
%     dros = cellfun(@(x) medfilt2(x,[5 5]),dros,'UniformOutput',0);
    tmp = cellfun(@(x) mean(x,2),dros,'UniformOutput',0);
    profile{1} = flipud(horzcat(tmp{:}));
    
    % smooth profile with a moving window of size 40 px
    y = smoothdata(profile{1},1,"movmean",40);
    x = ((1:size(y,1))-1).*tbl.pxSize(i);
    
    % plot intensity profile
    figure('Name',"s"+tbl.wells(i))
    hold on
    arrayfun(@(i)plot(x,y(:,i),'Color',line_color(i,:),'LineWidth',2),1:length(dros)) % 1:length(dros)
    hold off

    xlabel("Distance (mm)")
    ylabel("ROS (A.U.)")
    ylim([0 140])

    tick = get(gca,'YTick');
    ylimit = get(gca,'YLim');

    text(100,0.9*ylimit(end),"s"+tbl.wells(i)+", "+gapSize(i).*tbl.pxSize(i)+"µm")

    xline((0+[gapLeft(i)-1 gapRight(i)]).*tbl.pxSize(i),'--')

    set(gca,'Xlim',[0 2000],'XTick',0:500:2000,'XTickLabel',0:0.5:2,...
        'fontsize',16,'tickdir','out','XColor','k','YColor','k',...
        'LineWidth',1,'PlotBoxAspectRatio',[200 72 1])
    
    yyaxis right
    set(gca,'Ylim',[0 140],'YColor','k','YTickLabel',[])
    
    data = [data;repelem(round(gapSize(i).*tbl.pxSize(i)),length(x),1),x',y];

end

tbl_source = array2table(data,"VariableNames",["Gap size (um)" "Distance (um)" "T"+(1:7)]);
% writetable(tbl_source,"fig2a.xlsx")
