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

%% load kymomgraph
load('..\data\kymograph_s19.mat')

% ferroptosis
% data
ctrl.x = kymograph.front{1}.time;
ctrl.y = kymograph.front{1}.mode_dist;
ctrl.x = ctrl.x-ctrl.x(1);
ctrl.y = ctrl.y-ctrl.y(1);
% fitting line
p = polyfit(ctrl.x,ctrl.y,1);
ctrl.xf = [ctrl.x;20];
ctrl.yf = polyval(p,ctrl.xf);
ctrl.yf = ctrl.yf-ctrl.yf(1)/2;
ctrl.intersection = [9.4468 3.1063e+03];

% small molecule
sm.x = (0:0.1:25).*3600;
sm.y1 = sqrt(sm.x.*2.*200);
sm.y2 = sqrt(sm.x.*4.*200);
sm.y3 = sqrt(sm.x.*6.*200);
sm.intersection = [7.9714 3.3880e+03];

% protein
protein.x = sm.x;
protein.y1 = sqrt(protein.x.*2.*4.484);
protein.y2 = sqrt(protein.x.*4.*10);
protein.y3 = sqrt(protein.x.*6.*10);
protein.intersection = [14.5012 1.0218e+03];

xlimit = [0 20];
ylimit = [0 6000];

%% plot theoretical displacement
figure
hold on
%%% box
plot(xlimit,[ylimit(2) ylimit(2)],'k','LineWidth',1)
plot([xlimit(2) xlimit(2)],ylimit,'k','LineWidth',1)
plot(xlimit,[ylimit(1) ylimit(1)],'k','LineWidth',1)
plot([xlimit(1) xlimit(1)],ylimit,'k','LineWidth',1)
%%% line
p1 = plot(sm.x./3600,sm.y1,'k--','LineWidth',2);
p2 = plot(protein.x./3600,protein.y1,'k--','LineWidth',2);
p4 = plot(ctrl.xf,ctrl.yf,'r-','LineWidth',2);
p5 = plot(ctrl.x,ctrl.y,'ko','LineWidth',1,'MarkerFaceColor','w');
hold off

xlim(xlimit)
ylim(ylimit)
ylabel('Distance (mm)')
xlabel('Time (hr)')

lim = axis(gca);
set(gca,'YTickLabel',0:floor(lim(4)/1000),...
        'fontsize',16,'tickdir','out','XColor','k','YColor','k',...
        'LineWidth',1) % ,'PlotBoxAspectRatio',[16*(201+12)-12 ceil(5000/1.266) 1]

%%% annotation
pbaRatio = get(gca,'PlotBoxAspectRatio');
slope = p(1)/(diff(ylimit)/diff(xlimit))*(pbaRatio(2)/pbaRatio(1));
alpha = atand(slope);
text(ctrl.intersection(1)+1,ctrl.intersection(2),"Intercellular ferroptosis",...
    "HorizontalAlignment","center",'VerticalAlignment','top','Rotation',alpha,'FontSize',14)

slope = 200/sqrt(2*200.*sm.intersection(1)/3600)/(diff(ylimit)/diff(xlimit))*(pbaRatio(2)/pbaRatio(1));
alpha = atand(slope);
text(sm.intersection(1),sm.intersection(2),"Ca^{2+} diffusion",...
    "HorizontalAlignment","center",'VerticalAlignment','bottom','Rotation',alpha,'FontSize',14)

slope = 10/sqrt(2*10.*protein.intersection(1)/3600)/(diff(ylimit)/diff(xlimit))*(pbaRatio(2)/pbaRatio(1));
alpha = atand(slope);
text(protein.intersection(1),protein.intersection(2),"Protein diffusion",...
    "HorizontalAlignment","center",'VerticalAlignment','bottom','Rotation',alpha,'FontSize',14)

%%% add gradient background
lim = axis(gca);
xdata = [lim(1) lim(2)-0.05 lim(2)-0.05 lim(1)];
ydata = [lim(3) lim(3) lim(4)-0.05/21*6000 lim(4)-0.05/21*6000];
cdata(1,1,:) = [255 255 255]./255;
cdata(1,2,:) = [248 239 215]./255;
cdata(1,3,:) = [248 239 215]./255;
cdata(1,4,:) = [255 255 255]./255;
pp =patch(xdata,ydata,'k','Parent',gca);
set(pp,'CData',cdata, ...
    'FaceColor','interp', ...
    'EdgeColor','none');
uistack(pp,'bottom')
