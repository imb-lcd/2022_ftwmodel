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

%% load measurement data
tbl = readtable("..\data\Count.xlsx");
tbl = tbl(:,[1 2]);
tbl.Properties.VariableNames = ["gapSize" "isPassing"];

%% fit the data to a binary logistic model
x = tbl.gapSize;
y = tbl.isPassing;

% Fit model
mdl = fitglm(x, y, "Distribution", "binomial");
xnew = linspace(0,400,1000)'; % test data
ynew = predict(mdl, xnew);

ci = coefCI(mdl)

% sampling
isPlot = ismember(x,[70  46  75  77 111 111  54  80  39 101 187 193 139 160 150,...
    160 214 173 214 192 183 240 368 269 310 380 353 233 207 185]);

close all

figure('position',[680   250   560   560])

yyaxis left
scatter([x(isPlot);35;118;156;224],[y(isPlot);1;1;0;0],96,'ko','MarkerFaceAlpha',0.4);
yline([0 1],'--','Color',[0.6 0.6 0.6], 'LineWidth',2)
ylim([-0.1 1.1])
yticks([0 1])
yticklabels(["Halted" "Continued"])
set(gca,'YColor','k')

yyaxis right 
plot(xnew, ynew, 'r', 'LineWidth',2);
ylim([-0.1 1.1])
set(gca,'YColor','r')
set(gca,'fontsize',16.275,'TickDir','out','LineWidth',1.6275)
xlabel("Gap width (µm)",'FontSize',19.53)
ylabel("Propagating probability",'FontSize',19.53)
box off

% tbl = table([x(isPlot);35;118;156;224],[y(isPlot);1;1;0;0],'VariableNames',["Gap width (um)","isContinued"]);
% writetable(tbl,"../../source_data/fig2.xlsx","Sheet","b")

