%% initialize environment
clc; clear; close all

pos = "s19"; % position to be processed

% set imaging setting
frame_interval = 1; % hr
scale = 1.266; % um/px

%% for each image, estimate speed based on the kymograph
for i = 1

    % load mask of initiation sites
    load("..\data\setting_"+pos(i)+".mat")
    
   % load images for speed estimation
    list = string(ls("..\img\*10.tif"));
    list = list(contains(list,pos(i)) & contains(list,["cy5" "ros"])); % 2 signals are allowed
    img = cell(size(list));
    for j = 1:length(list)
        info = imfinfo("..\img\"+list(j));
        img{j} = arrayfun(@(i)imread("..\img\"+list(j), i, 'Info', info),...
            1:length(info),'UniformOutput',0);
    end
    
    % set directions for speed measurement
    % determine the directions based on the location of an initiation
    K = kron([3 4 5;2 9 6;1 8 7],ones(ceil(size(initiation.mask)./3)));
    K = K(1:size(initiation.mask,1),1:size(initiation.mask,2));
    tbl = tabulate(K(initiation.mask));
    [~,index] = max(tbl(:,2));
    if index<9
        rotation_angle = mod((0:5:90)+(45.*(index-1)),360);
    else
        rotation_angle = 0:5:359;
    end
    % or directly specify a direction by an angle
    rotation_angle = 40; 
    
    % calculate the initiation center after rotation
    tmp = arrayfun(@(a)regionprops('table',imrotate(initiation.mask,a),{'Centroid' 'area'}),...
        rotation_angle,'UniformOutput',0);   
    center = nan(length(tmp),2); % [x,y]
    for j = 1:length(tmp)
        [~,index] = max(tmp{j}.Area);
        center(j,:) = fix(tmp{j}.Centroid(index,:));
    end
    
    % identify the furthest index can be used
    tmp = arrayfun(@(a)imrotate(ones(size(initiation.mask)),a),rotation_angle,'UniformOutput',0);
    tmp = arrayfun(@(i)tmp{i}(1:center(i,2),center(i,1)-100:center(i,1)+100),1:length(tmp),'UniformOutput',0);
    index_strt = cellfun(@(a)find(sum(a==1,2)==201,1),tmp);
    
    % kymograph, to keep speed measurement for each direction
    [C1,C2]=meshgrid(1:length(list),1:length(rotation_angle));
    kymograph = table(list(C1(:)),rotation_angle(C2(:))',C1(:),C2(:),...
        'VariableNames',["source_img","angle","index_img","index_angle"]);
    
    % get kymographs for each directions
    for j = 1:height(kymograph)
        imgRotated = imrotate3(cat(3,img{kymograph.index_img(j)}{:}),...
            kymograph.angle(j),[0 0 1]);
        imgCropped = imgRotated(1:center(kymograph.index_angle(j),2),...
            center(kymograph.index_angle(j),1)-100:center(kymograph.index_angle(j),1)+100,:);
        maxR = max(imgCropped,[],2);
        kymograph.maxR{j} = rot90(maxR(index_strt(kymograph.index_angle(j)):end,:),-1);
    end

    % speed measurement
    th = 0.042; % a threshold to identify the propagation front
    for j = 1:height(kymograph)
        maxR = kymograph.maxR{j}; % kymograph

        % remove all-black images if necessary
        index_frame = (1:size(maxR,1))'; % record the remaining frame
        index_noImg = find(sum(maxR==0,2)==size(maxR,2));
        if ~isempty(index_noImg)
            maxR(index_noImg,:) = [];
            index_frame(index_noImg) = [];
        end

        % identify the propagation front
        X = imbinarize(imadjust(maxR,stretchlim(maxR)),th);
        Y = X .* (index_frame);  % Replace 1s by the column index
        Y(X == 0) = NaN;            % Let all but the former 1s be NaN
        % if there are at least 2 frames are above threshold, mark the
        % pixel as reached by waves
        for k = 1:size(Y,2)
            if nnz(~isnan(Y(:,k)))>2
                M = max(Y(:,k));
                X((M+1):end,k) = 1;
            end
        end
        [G,ID] = findgroups(min(Y, [], 1));
        G = arrayfun(@(x)find(G==x),ID,'UniformOutput',0);
        G = table(ID',G','VariableNames',{'frame' 'px'});

        % idenitfy the range of frames for speed estimation
        test = sum(X>0,2)./size(X,2).*100;
        index = find(normalize(test,'range')>0.15 & normalize(test,'range')<0.85);
        frame_strt = index_frame(index(find(diff(index,2)==0,1))); % identify the beginning of a consecutive period
        frame_stop = index_frame(index(find(diff(index,2)==0,1,'last')+2)); % identify the end of a consecutive period
        frame_strt=6;frame_stop=18; % or directly specify the frames for speed estimation

        % keep the wavefront within the period
        G = G(ismember(G.frame,frame_strt:frame_stop),:);
        G(cellfun('length',G.px)==0,:) = [];

        % convert frames and pixels to real time and distance
        G.time = (G.frame-1).*frame_interval;
        G.dist = cellfun(@(x)(x-1).*scale,G.px,'UniformOutput',0);
        
        % find the neighbors around the mode of distance for each frame
        nNbr = 11;
        [f,xi] = cellfun(@(a)ksdensity(a),G.dist,'UniformOutput',0);
        [~,I] = cellfun(@(a)max(a),f);
        x_mode = arrayfun(@(i)xi{i}(I(i)),1:length(I))';
        [~,I] = arrayfun(@(i)sort(abs(G.dist{i}-x_mode(i))),1:length(I),'UniformOutput',0);
        G.data_dist = arrayfun(@(i)G.dist{i}(I{i}(1:min([nNbr length(G.dist{i})]))),1:length(I),'UniformOutput',0)';
        G.mode_dist = x_mode;

        % fit a straight line based on the wavefront data
        p = cell(size(frame_strt,1),1);
        for k = 1:size(frame_strt,1)
            data_t = repelem(G.time,cellfun('length',G.data_dist));
            data_x = [G.data_dist{:}]';
            p{k} = polyfit(data_t,data_x,1);
        end

        % organize data
        kymograph.front{j} = G; 
        kymograph.frame_range{j} = [frame_strt frame_stop]; 
        kymograph.speed{j} = p;

    end
    
    % plot kymograph
    figure
    tiledlayout('flow','TileSpacing','compact','Padding','compact')
    order = reshape(1:height(kymograph),length(rotation_angle),[])';
    for j = order(:)'
        if ~isempty(kymograph.front{j}) % only plot when speed is measured
            nexttile
            hold on
            plotKymo(kymograph.maxR{j},frame_interval,scale)
            plot(kymograph.front{j}.time,kymograph.front{j}.mode_dist,'co')
            hold off

            sig = extractBetween(kymograph.source_img(j),"_","_");
            title(sprintf("%s,%d,%.2f",sig,kymograph.angle(j),kymograph.speed{j}{1}(1)/60))
            xticklabels("")
            yticklabels("")
        end
    end

end

