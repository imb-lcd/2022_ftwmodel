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

clearvars
fname = ["DFO_24_s7" "DFO_24_s14" "FAC_24_s1" "FAC_24_s5",...          %  1  2  3 4
    "LY294002_24_s19" "LY294002_24_s21",...                            %  5  6
    "Dasatinib_20211113_s24" "Dasatinib_20211117_s4",...               %  7  8
    "GKT_24_20210827_s9" "GKT_24_20210827_s14" "GKT_24_20211113_s7",...%  9 10 11
    "GKT_24_20211113_s8" "GKT_24_20211117_s11",...                     % 12 13
    "DFO_24_s11_c1" "DFO_24_s13_c1" "DFO_24_s23_c1",...                % 14 15 16
    "FAC_24_s2_c1" "FAC_24_s6_c1",...                                  % 17 18
    "LY294002_24_s4_c1" "LY294002_24_s8_c1" "LY294002_24_s9_c1",...    % 19 20 21
    "LY294002_24_s16_c1" "LY294002_24_s17_c1" "LY294002_24_s18_c1",... % 22 23 24
    "Dasatinib_24_20211117_s9_c1" "Dasatinib_24_20211117_s15_c1",...   % 25 26
    "Dasatinib_24_20211113_s20_c1",...                                 % 27
    ]+"_ff10";

frame1_treat = [12 12 16 16,...
    12 12,...
    13 12,...
    12 12 14,...
    14 12,...
    12 12 12,...
    16 16,...
    12 12 12,...
    13 13 13,...
    12 12,...
    12,...
    ];
%%
for well = [2 4 6 10 8]
    oname = "..\data\setting_"+fname(well)+".mat";
    load(oname)
    frame_strt = speed.cy5.frame_strt;
    frame_stop = speed.cy5.frame_stop;

    kymo = "kymograph_cy5_"+fname(well)+".mat";
    load("..\data\"+kymo)

    frame_interval = 1;
    scale = 1.266;
    if ismember(well,[6 8])
        scale = 1.260;
    end
    offset = 0; % time-shift the kymograph

    maxR = kymograph{1,1};

    % x- and y-limit
    xlimit = [0-0.5 size(maxR,1)-1+0.5].*frame_interval;
    ylimit = [0 size(maxR,2)-1].*scale;

    % background
    len_px = max(size(maxR));
    B = imresize(maxR,[len_px len_px]);
    H = fspecial('average',11);
    MotionBlur = imfilter(B,H,'replicate');
    MotionBlur = imadjust(uint16(MotionBlur),stretchlim(uint16(MotionBlur),[0.01 0.99]));

    % plot
    figure

    % kymograph
    hold on
    imagesc([0-0.5 size(maxR,1)-1+0.5].*frame_interval,...
        [size(maxR,2)-1 0].*scale,MotionBlur')
    hold off
    colormap gray
    xlim(xlimit)
    ylim(ylimit)

    set(gca,'Xlim',xlimit,'YLim',ylimit,'YTickLabel',0:floor(ylimit(2)/1000),...
        'fontsize',16,'tickdir','out','XColor','k','YColor','k',...
        'LineWidth',1,'Layer','top','YTick',0:1000:floor(ylimit(2)/1000)*1000)
    xlabel('Time (hr)')
    ylabel('Distance (mm)')
    hold off
    axis square

    % add top and left edge
    hold on
    plot(xlimit,[ylimit(2) ylimit(2)],'k','LineWidth',1)
    plot([xlimit(2) xlimit(2)],ylimit,'k','LineWidth',1)

    % label treatment time
    plot(frame1_treat(well)-0.5-offset,0.9*ylimit(2),'wv','MarkerFaceColor','w')
    hold off

    % label condition
    text(0.5,0.9*ylimit(2),strrep(fname(well),'_','\_'),'color','w','Rotation',-90)

    % fitting line
    displacement = 0;
    hold on
    for i = 1:size(frame_strt)
        p = kymograph{1,5}{i};
        slope = p(1)/(diff(ylimit)/diff(xlimit)); % should match the xlimit and ylimit of figure
        ang = atand(slope);

        xrange = (frame_strt(i,1):frame_stop(i,1))-1;
        F = polyval(p,xrange);
        plot(xrange-offset,F+slope*displacement,...
            'LineWidth',2,'color',[1 1 1 0.7])

        v1 = [frame_strt(i,1)-offset-1,F(1)]; %(x,y)
        v2 = [frame_stop(i,1)-offset-1,F(end)]; %(x,y)

    end
    hold off
    close all

end