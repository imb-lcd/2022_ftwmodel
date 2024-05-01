%% initialize environment
clc; clear; close all

isSave = 1; % 1: save the result; 0: don't save the result
% condition = "RPE_1_Erastin";
well = "11";
frame_range = 1:5;

%% prepare landscape
% load outlines of waves
load("..\data\mask_s"+well+".mat","MASK")
mask = MASK(frame_range);

% create a landscape
landscape = uint8(cat(3,mask{:}));
[a,Z]=max(landscape,[],3);
Z(a==0) = size(landscape,3)+1;

% calculate gradient
sz = 5; % every sz grid points get a vector. Smaller value will take longer time.
[Xq,Yq] = meshgrid(1:sz:size(Z,1));
Vq = Z(1:sz:size(Z,1),1:sz:size(Z,2));
[U,V] = gradient(Vq); % generate gradient based on the all frames of interest

% solve the conflicts of vectors' direction (should ignore length average?)
UG = imgaussfilt(U,3);
VG = imgaussfilt(V,3);
arrow_angle = atan2(VG,UG);

% determine if the arrows are valid for entropy calculation
isValid = true(size(Xq));
isValid(Vq==max(Z(:))) = 0; % within the frame range (frame 1 to 5)
isValid(sqrt(UG.^2+VG.^2)==0) = 0; % length of arrow > 0

%% set parameters of bounding box (initiations and size)
if exist("..\data\center_s"+well+".mat","file")
    load("..\data\center_s"+well+".mat")
else
    error("No initiation is found.")
end
pgon = arrayfun(@(i)polyshape(center{i}(1)+[-1 -1 1 1 -1].*ROI_width./2,...
    center{i}(2)+[0 -1 -1 0 0].*ROI_length),1:length(center),'UniformOutput',0);

%% collect information of valid arrows for each initiation
% set the angles of view
alpha = (0:20:359)'; % the angle of view, up=0, right=90, down=180, left=270   

% identify the position of valid arrows
valid_grid = cell(size(alpha));
for i = 1:length(alpha)

    % rotate the bounding box
    polyout = cellfun(@(p,c)rotate(p,alpha(i),c'),pgon,center,'UniformOutput',0);

    % identify the arrow within the rotated bounding box
    in = cellfun(@(p)reshape(inpolygon(Xq(:),Yq(:),p.Vertices(:,1),p.Vertices(:,2)),size(Xq,1),[]),...
        polyout,'UniformOutput',0);
    [a,BC] = max(uint8(reshape([in{:}],size(in{1},1),size(in{1},2),length(in))),[],3);
    BC(a==0) = 0;
    valid_grid{i} = isValid.*BC;
    
    % check if the arrows obtained from at least 3 frames
    index = find(arrayfun(@(j)length(unique(Vq(valid_grid{i}==j))),1:length(center))<3);
    if ~isempty(index)
        valid_grid{i}(ismember(valid_grid{i},index)) = 0;
    end

    % check if the number of valid arrows is greater than 20
    index = find(arrayfun(@(j)nnz(valid_grid{i}==j),1:length(center))<=20);
    if ~isempty(index)
        valid_grid{i}(ismember(valid_grid{i},index)) = 0;
    end
    
end

% collect valid arrows
data_arrow = [];
for i = 1:length(valid_grid)
    index = arrayfun(@(j)find(valid_grid{i}==j),1:length(center),'UniformOutput',0);
    
    % get rid of 20% of arrows, which are deviated from the angle of view
    % if the number of remaining arrows is less than 100, keep all arrows; 
    % otherwise, use 80% of the arrows
    tmp = cellfun(@(c)angdiff(arrow_angle(c),repelem(deg2rad(alpha(i)-90),length(c),1)),...
        index,'UniformOutput',0);
    [~,II] = cellfun(@(a)sort(abs(a),'ascend'),tmp,'UniformOutput',0);
    for j = 1:length(II)
        if length(II{j})>100
            II{j} = II{j}(1:fix(0.8*length(II{j})));
        end
    end
    index = cellfun(@(a,b)a(b),index,II,'UniformOutput',0);
    
    % collect valid arrows
    nArrow = cellfun('length',index);
    index = vertcat(index{:});
    tmp = [repelem(alpha(i),length(index),1),Xq(index),Yq(index),...
        UG(index),VG(index),arrow_angle(index),repelem(1:length(center),nArrow)'];
    data_arrow = [data_arrow;tmp];
end
data_arrow = array2table(data_arrow,'VariableNames',["angleBox" "x" "y" "u" "v" "angleArrow" "centerID"]);

%% plot polarhistogram and calculate entropy
binedges = linspace(0,2*pi,31);
H = nan(length(alpha),length(center));
for i = 1:length(center)
    figure('Name',"Center "+i,'Position',[50 50 960 480])
    tiledlayout(3,6,'TileSpacing','compact','Padding','compact')
    for j = 1:length(alpha)
        a = data_arrow.angleArrow(data_arrow.angleBox==alpha(j)&data_arrow.centerID==i);
        if ~isempty(a)
            nexttile(j)
            h = polarhistogram(-1.*a,'BinEdges',binedges,'Normalization','probability');
            p = h.Values;
            p(p==0)=[]; % remove the zero probability
            
            H(j,i) = -sum(p.*log2(p)); % calculate the entropy
            
            title("S="+sprintf("%.3f",H(j,i))+", "+alpha(j)+"^\circ")
            thetaticks([0 45 90 135 180 225 270 315])
            thetaticklabels({'90','45','0','315','270','225','180','135'})
        end
    end
end
