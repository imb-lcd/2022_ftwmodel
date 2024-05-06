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

%% prepare the filename to be analyzed
wells = 1:60;
wells = reshape(wells,12,[])';
wells(2,:) = fliplr(wells(2,:));
wells(4,:) = fliplr(wells(4,:));
wells = wells(:);

conc = [160./2.^(0:10) 0]; % 160/2^11 is 0
level = 1:length(conc);
conc = repmat(conc,5,1);
conc = conc(:);

level = fliplr(level);
level = repmat(level,5,1);
level = level(:);

fname = "kymograph_cy5_DFO_96_s"+(1:60)+"_c1_ff10.mat";
frame1_treat = 4.*ones(size(fname));

frame_interval = 1; % hr
scale = 1.260; % um

%% organize the measured speed
speeds = nan(length(wells),1);
angle = speeds;
pct = speeds;
frame_strt = speeds;
frame_stop = speeds;
for i = 1:length(wells)
    load("..\data\"+fname(wells(i)))
    speeds(i) = abs(kymograph{5}{1}(1))/3600;
    
    load("..\data\"+strrep(fname(wells(i)),"kymograph_cy5","setting"))
    angle(i) = speed.rotation_angle;
    pct(i) = speed.cy5.pct;
    frame_strt(i) = speed.cy5.frame_strt;
    frame_stop(i) = speed.cy5.frame_stop;

end

result = table(wells,conc,level,speeds,angle,pct,frame_strt,frame_stop);
result(40,:) = [];

%% fit a curve for the dose-speed relationship
y = result.speeds.*60;%(result.isoutlier==0);
well = result.wells;%(result.isoutlier==0);

% get the real concentration of the lowest level
x = result.conc;

% set speed of 160uM as 0
y(x==160) = 0;
isUsed = true(size(x));
x = x(isUsed);
y = y(isUsed);
well = well(isUsed);

% add 5 zeros to 320/640uM
nZero = 0;
x = [x;repmat(160.*ones(nZero,1),5,1)];
y = [y;repmat(zeros(nZero,1),5,1)];
well = [well;repmat(100.*(1:nZero)',5,1)];

% model of a dose-response curve
F = @(z,xdata) z(1)+(z(2)-z(1)).*z(3)./(z(3)+xdata);
z0 = [0  6   80];
lb = [0  0    0];
ub = [6  9  Inf];
[z,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(F,z0,x,y,lb,ub);%
ci = nlparci(z,residual,'jacobian',jacobian);

fx = 0:min(abs(diff(unique(x)))):max(x);

[ypred,delta] = nlpredci(F,fx,z,-1.*residual,'Jacobian',jacobian,'alpha',0.05);
lower = ypred - delta;
upper = ypred + delta;
lower(lower<0)=0;

% take equally spaced points on log10 axis
xq = 10.^(linspace(log10(min(unique(x(x>0)))),log10(max(unique(x(x>0)))),101));
lower = interp1(fx,lower,xq,'spline');
upper = interp1(fx,upper,xq,'spline');
ypred = interp1(fx,ypred,xq,'spline');
fx = xq;

%% plot dose-response curve
close all
figure('position',[1020 500 420 420])
t = tiledlayout(1,34,'TileSpacing','tight');
bgAx = axes(t,'XTick',[],'YTick',[],'Box','off');
axis square
bgAx.Layout.TileSpan = [1 34];

ax1 = axes(t);
ax1.Layout.Tile = 1;
ax1.Layout.TileSpan = [1 5];
p = plot(x(x==0),y(x==0),'o','color',[59, 73, 146]/255,'MarkerFaceColor','none');
row = dataTipTextRow('Well',well(x==0));
p.DataTipTemplate.DataTipRows(end+1) = row;
set(ax1,'fontsize',16,'YLim',[0 6.5],'LineWidth',1,'XTick',0,'XTickLabel',0,'Box','off','XLim',[-0.1 0.1])

ax2 = axes(t);
ax2.Layout.Tile = 6;
ax2.Layout.TileSpan = [1 34-ax2.Layout.Tile+1];
ax2.YAxis.Visible = 'off';
set(ax2,'Box','off','XScale','log','XLim',[0.05 200],'YLim',[0 6.5],'FontSize',16,'LineWidth',1,'XTick',10.^(-1:2),'XTickLabel',10.^(-1:2))
hold on
patch([fx,fliplr(fx)],[lower,fliplr(upper)],[.6 .6 .6],'facecolor',[.9 .9 .9],'edgecolor','none')
plot(fx,ypred,'k--','LineWidth',1)
plot(fx,[lower;upper],'LineWidth',1,'color',[.6 .6 .6])
p = plot(x(x>0),y(x>0),'o','color',[59, 73, 146]/255,'MarkerFaceColor','none');
row = dataTipTextRow('Well',well(x>0));
p.DataTipTemplate.DataTipRows(end+1) = row;
hold off

% Link the axes
linkaxes([ax1 ax2], 'y')
set(bgAx,'Color','w','XColor','none')
xlabel(t,"DFO (µM)",'fontsize',16)
ylabel(t,"Speed (µm/min)",'fontsize',16)

text(0.1,3,'$$y=y_0+(y_M-y_0)\frac{K}{K+x}$$',"Interpreter","latex")
text(0.1,2.5,sprintf('$$y_0=%.3f\\pm%.3f$$',z(1),z(1)-ci(1,1)),"Interpreter","latex")
text(0.1,2,sprintf('$$y_M=%.3f\\pm%.3f$$',z(2),z(2)-ci(2,1)),"Interpreter","latex")
text(0.1,1.5,sprintf('$$K=%.3f\\pm%.3f$$',z(3),z(3)-ci(3,1)),"Interpreter","latex")

