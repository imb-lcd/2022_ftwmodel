%% initialize environment
clc
clear
close all

pos = "s"+[41,1,5,12,15,16,29,36]; % well to be processed

%% for each well, identify one initiation
for i = 1:length(pos)

    % set the filename of the output
    oname = "..\data\setting_"+pos(i)+".mat";

    % if the output file existed, load it; otherwise use default setting
    if exist(oname,'file')
        load(oname)
        isSave = 0;
    else
        % default setting
        % frame_index: frame of initiation
        % signal: which signal is going to be used for identify initiation
        % range_area: filter out the region whose area is outside the range
        % r_dilation/r_erosion: connect dead cells
        % th_mask: the area of connected regions smaller than th_mask will
        %          be removed.
        % boundary_index: the index of the region where the initiation
        %                 locates
        initiation = struct('frame_index',3,'signal',"cy5",'th_intensity',0.9,...
            'range_area',[20 1000],'r_dilation',15,'r_erosion',10, ...
            'th_mask',5000,'boundary_index',1);
        isSave = 1;
    end
    
    % load image
    imgname = "..\img\"+pos(i)+"_"+initiation.signal+"_ff1_10.tif";
    img = imread(imgname, initiation.frame_index, 'Info', imfinfo(imgname));
    
    % binarize image
    im_adjusted = imadjust(img,stretchlim(img));    
    mask = imbinarize(im_adjusted,initiation.th_intensity);

    % filter out the regions with abnormal area
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

    figure('WindowState','maximized')
    tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
    nexttile
    imshow(L,[0 max(L(:))],'Border','tight');
    title('Label matrix')    
    nexttile
    imshow(B)
    title('Image overlaid with the label matrix')
    
    % mark the initiation
    maskDS = L==initiation.boundary_index;    
    hold on
    stats = regionprops(maskDS,'Centroid');
    plot(stats(1).Centroid(1),stats(1).Centroid(2),'rx','MarkerSize',20)

    % save mask into setting
    initiation.mask = maskDS;
    if isSave
        save(oname,"initiation")
        exportgraphics(gcf,"..\img\initiation_"+pos(i)+".jpg","Resolution",300)
    end
    
    close all
end

%%
% % pos = "s"+[41,1,5,12,15,16,29,36]; % well to process
% clear
% load('..\data\setting_s36.mat')
% new = initiation;
% load('..\data\v1\setting_s36.mat')
% old = initiation;
% imshowpair(old.mask,new.mask)