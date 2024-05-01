%% initialize environment
clc
clear
close all

pos = "s"+[19,4,8,9,10,21,22,23]; % well to be processed
angle = [40 75 210 30 355 90 85 40]; % the direction to plot kymograph
angle = [40 85 210 0 0 75 90 40]; % the direction to plot kymograph
% angle = [90 80 215 90 345 70 70 20];

% set the parameter for visualizing kymograph
frame_interval = 1;
scale = 1.266;
frame_toshow = 4:19; % frame range to show

xlimit = ([0 frame_toshow(end)-frame_toshow(1)]+[-0.5 0.5]).*frame_interval; 
ylimit = [0 5000];
pbaRatio = [length(frame_toshow)*(201+12)-12 ceil(ylimit(2)/scale) 1];

color_signal = [255 255 255;254 212 57]./255; % white/yellow

%% plot kymograph
speed = [];
for i = 1

    % load setting for setting and kymograph
    load("..\data\setting_"+pos(i)+".mat");
    load("..\data\kymograph_"+pos(i)+".mat");

    % get kymograph
    tmp = kymograph(ismember(kymograph.angle,angle(i)),:);

    % get speed
    speed = [speed;tmp.speed{1}{1}(1)/60];

    % retrieve kymograph
    maxR = tmp.maxR{1};
    frame_range = 1:size(maxR,1);
    
    % calculate the offset of x- and y-axis (spatial and temporal axis)
    offset_y = 0;%(size(maxR,2)-ceil(ylimit(2)/scale)).*scale;
    offset_x = (min(frame_toshow)-1).*frame_interval;

    maxR = maxR(frame_toshow,1:ceil(ylimit(2)/scale));
    
    % adjust kymograph of well 19 for a dimer frame (frame 16) by
        % scaling the intensity to the mean of neighboring frames (frame 15
        % and 17).
    if i == 1
        adj_pct = 50;
        maxR_adjusted = double(maxR);        
        maxR_adjusted(13,:) = maxR_adjusted(13,:).*...
            prctile(mean(maxR_adjusted([12,14],1:3000))./maxR_adjusted(13,1:3000),adj_pct);
        maxR = maxR_adjusted;
    end

    % plotting
    figure
    hold on
    plotKymo(maxR,frame_interval,scale)
    
    % fitting line - cy5
    F = polyval(tmp.speed{1}{1},(frame_range-1).*frame_interval)-offset_y;
    plot((frame_range-1).*frame_interval-offset_x,F,'LineWidth',2,'color',[color_signal(1,:) 0.7]) % cy5

    % annotate the measured speed - cy5
    x = (frame_range-1).*frame_interval-offset_x;
    y = F;
    [xi,yi] = polyxpoly(xlimit,[mean(ylimit) mean(ylimit)],x,y);
    slope = tmp.speed{1}{1}(1)/(diff(ylimit)/diff(xlimit))*(pbaRatio(2)/pbaRatio(1));
    alpha = atand(slope);
    text(0,4500,pos(i),'color','w')
    text(xi-1,yi,round(tmp.speed{1}{1}(1)/60,2)+"µm/min",...
        'color','w','HorizontalAlignment','center','Rotation',alpha,...
        'VerticalAlignment','bottom','FontSize',16)

    % set axis range and labels
    set(gca,'Xlim',xlimit,'YLim',ylimit,'YTickLabel',0:5,...
        'fontsize',16,'tickdir','out','XColor','k','YColor','k',...
        'LineWidth',1,'PlotBoxAspectRatio',pbaRatio,'Layer','top')
    
    % add top and left edge
    plot(xlimit,[ylimit(2) ylimit(2)],'k','LineWidth',1)
    plot([xlimit(2) xlimit(2)],ylimit,'k','LineWidth',1)

    xlabel('Time (h)')
    ylabel('Distance (mm)')

    hold off

end
