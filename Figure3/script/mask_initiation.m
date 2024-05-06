%% initialize environment
clc
clear
close all

% available drug: DFO, FAC, Dasatinib, LY294002, GKT
% 24/96 wells
drug = "GKT";
plate = 24;

% find related image filename
list = string(ls("..\img\"+drug+"_"+plate+"_*.tif"));
expression = ['(?<drug>\w+)_(?<plate>\d+)_s(?<well>\d+)_c(?<ch>\d+)_ff1_10.tif'];
tokenNames = regexp(list,expression,'names');
pos = struct2table([tokenNames{:}]');

%% for each well, identify one initiation
for i = 1:height(pos)

    % set the filename of the output
    if exist("..\data\"+pos.drug(i)+"_"+pos.plate(i)+"\","dir")==0
        mkdir("..\data\"+pos.drug(i)+"_"+pos.plate(i))
    end
    oname = "..\data\"+pos.drug(i)+"_"+pos.plate(i)+"\setting_s"+pos.well(i)+".mat";

    % if the output file existed, load it; otherwise use default setting
    if exist(oname,'file')==2
        load(oname)
        isSave = 1;
    else
        % default setting
        % frame_index: frame of initiation
        % signal:
        % range_area: filter out the region whose area is outside the range
        % r_dilation/r_erosion: connect dead cells
        % th_mask: 
        % boundary_index:

        initiation = struct('frame_index',2,'signal',"cy5",'th_intensity',0.9,...
            'range_area',[20 5000],'r_dilation',15,'r_erosion',10, ...
            'th_mask',5000,'boundary_index',1);
        clearvars tmp
        isSave = 1;
    end
    
    % load image
    imgname = "..\img\"+list(i);
    img = imread(imgname, initiation.frame_index, 'Info', imfinfo(imgname));
    
    % binarize image
    im_adjusted = imadjust(img,stretchlim(img));    
    mask = imbinarize(im_adjusted,initiation.th_intensity);
    mask = bwareafilt(mask,initiation.range_area);
    % imshow(BW)

    % dilate the threholding result to connect dead cells
    maskD = imfill(bwdist(mask)<initiation.r_dilation,'holes');    
    maskDS = bwdist(~maskD)>initiation.r_erosion;

    % remove small objects from the mask
    maskDS = bwareaopen(maskDS,initiation.th_mask);  

    % examine the boundary
    L = bwlabel(maskDS); 
    B = labeloverlay(im_adjusted,bwlabel(maskDS));
    
    % check initiation
    figure("Name",pos.drug(i)+"_"+pos.plate(i)+"_s"+pos.well(i),'WindowState','maximized')
    tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
    nexttile
    imshow(L,[0 max(L(:))],'Border','tight');
    title('Label matrix')    
    nexttile
    imshow(B)
    title('Image overlaid with the label matrix')
    
    % mark the initiation
    maskDS = L==initiation.boundary_index;
    initiation.mask = maskDS;
    hold on
    stats = regionprops(maskDS,'Centroid');
    plot(stats(1).Centroid(1),stats(1).Centroid(2),'rx','MarkerSize',20)

    % save mask into setting    
    if isSave
        save("..\data\"+oname,"initiation")
        exportgraphics(gcf,"..\img\initiation_check\initiation_"+...
            pos.drug(i)+"_"+pos.plate(i)+"_s"+pos.well(i)+".jpg",...
            "Resolution",300)
    end
    
    close all
end

%%
% pos = "s"+[19,4,8,9,10,21,22,23]; % well to process
% clear
% load('..\data\setting_s23.mat')
% new = initiation;
% load('..\data\v1\setting_s23.mat')
% old = initiation;
% imshowpair(old.mask,new.mask)