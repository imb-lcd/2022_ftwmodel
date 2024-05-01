%% initialize environment
clc; clear; close all

isSave = 0; % 1: save the result; 0: don't save the result
% condition = "RPE_1_Erastin";
well = "11";
frame_range = 1:5;
pixelsize = 1.266;

%% prepare landscape
% load outlines of waves
load("..\data\mask_s"+well+".mat","MASK")
mask = MASK(frame_range);

% create a landscape
landscape = uint8(cat(3,mask{:}));
[a,Z]=max(landscape,[],3);
Z(a==0) = size(landscape,3)+1;

% calculate gradient based on the landscape
sz = 20; % every sz grid points get a vector. Smaller value will take longer time.
[Xq,Yq] = meshgrid(1:sz:size(Z,1));
Vq = Z(1:sz:size(Z,1),1:sz:size(Z,2));
[U,V] = gradient(Vq); % generate gradient based on the all frames of interest

% solve the conflicts of vectors' direction (should ignore length average?)
UG = imgaussfilt(U,1);
VG = imgaussfilt(V,1);

% set the same length for all vectors
L = sqrt(UG.^2+VG.^2);
UG = sz.*UG./L;
VG = sz.*VG./L;

%% set parameters of bounding box (initiations and size)
if exist("..\data\center_s"+well+".mat","file")
    load("..\data\center_s"+well+".mat")
else
    ROI_width = 75;
    ROI_length = 275; 
    center{1} = [399;1325];
    center{2} = [2425;1207];
    center{3} = [2286;2464];
end

%% global view of vector field
% plot laandscape
figure;imshow(Z-0.5,[min(Z(:))-1 max(Z(:))],'Border','tight')
colormap(gray(max(Z(:))-min(Z(:))+1))
% copygraphics(gcf)

% show colorbar for different frames
% color_frame = gray(6);
% color_frame(end,:) = [];
% colormap(color_frame)
% caxis([0 5])
% colorbar('ticks',(1:5)-0.5,'TickLabels',1:5)
% colorbar('ticks',unique(Z(:))-0.5,'TickLabels',unique(Z(:)))

% plot arrows
hold on
color_vector= colorcet('R4','N',max(Z(:))-1,'reverse',1);
arrayfun(@(i)quiver(Xq(Vq==i),Yq(Vq==i),UG(Vq==i),VG(Vq==i),0,...
    'color',color_vector(i,:)),1:length(mask));

% mark initiation and bounding boxes of entropy calculation
for i = 1:length(center)
    pgon = rotate(nsidedpoly(8,'Center',center{i}','Radius',ROI_length),45/2,center{i}');
    x = [reshape(pgon.Vertices(:,1),4,[]),nan(4,1)]';
    y = [reshape(pgon.Vertices(:,2),4,[]),nan(4,1)]';
    plot(x(:),y(:),'w--','LineWidth',2)
    text(center{i}(1),center{i}(2),string(i),'fontweight','bold',...
        'FontSize',14,'Color','m','HorizontalAlignment','center','VerticalAlignment','middle')
end

%% local view of vector field (around initiations)
blankspace = ceil(ROI_length*(1/0.85-1)/10)*10;
xlimit = cellfun(@(c)[max([c(1)-ROI_length-blankspace 1]) min([c(1)+ROI_length+blankspace size(Z,2)])],center,'UniformOutput',0);
ylimit = cellfun(@(c)[max([c(2)-ROI_length-blankspace 1]) min([c(2)+ROI_length+blankspace size(Z,1)])],center,'UniformOutput',0);
side_length = cellfun(@(c)min([c(1)-1 size(Z,2)-c(1) c(2)-1 size(Z,1)-c(2)]),center);

[r,c] = find([vertcat(xlimit{:}),vertcat(ylimit{:})]== ...
    repmat([1 size(Z,2) 1 size(Z,1)],length(center),1));

% if the default axis limit touches the edge, using the minimal distance as
% sidelength to determine the distance
if ~isempty(r)
    for i = 1:length(r)
        xlimit{r(i)} = [max([center{r(i)}(1)-side_length(r(i)) 1]) min([center{r(i)}(1)+side_length(r(i)) size(Z,2)])];
        ylimit{r(i)} = [max([center{r(i)}(2)-side_length(r(i)) 1]) min([center{r(i)}(2)+side_length(r(i)) size(Z,1)])];
    end
end

for j = 1:length(xlimit)
    % plot laandscape
    figure;imshow(Z-0.5,[min(Z(:))-1 max(Z(:))],'Border','tight')
    colormap(gray(max(Z(:))-min(Z(:))+1))

    % plot arrows
    hold on
    color_vector= colorcet('R4','N',max(Z(:))-1,'reverse',1);
    arrayfun(@(i)quiver(Xq(Vq==i),Yq(Vq==i),UG(Vq==i),VG(Vq==i),0,...
        'color',color_vector(i,:)),1:length(mask));

    % mark initiation and bounding boxes of entropy calculation
    for i = 1:length(center)
        pgon = rotate(nsidedpoly(8,'Center',center{i}','Radius',ROI_length),45/2,center{i}');
        x = [reshape(pgon.Vertices(:,1),4,[]),nan(4,1)]';
        y = [reshape(pgon.Vertices(:,2),4,[]),nan(4,1)]';
        plot(x(:),y(:),'w--','LineWidth',2)
        text(center{i}(1),center{i}(2),string(i),'fontweight','bold',...
            'FontSize',14,'Color','m','HorizontalAlignment','center','VerticalAlignment','middle')
    end

    % set axit limit
    xlim(xlimit{j})
    ylim(ylimit{j})

    % add scale bar (200 um)
    hold on    
    pgon = polyshape(xlimit{j}(1)+10+[0 200/pixelsize 200/pixelsize 0 0],...
        ylimit{j}(1)+10+[0 0 20 20 0]);
    p = plot(pgon);
    p.EdgeColor='none';
    p.FaceAlpha=1;
    p.FaceColor='k';
    hold off
    text(xlimit{j}(1)+10+100/pixelsize,ylimit{j}(1)+10+30,"200 um",'HorizontalAlignment','center')
    
end
