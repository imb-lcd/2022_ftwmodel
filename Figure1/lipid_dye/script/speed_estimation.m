%% initialize environment
clc
clear
close all

pos = "s"+[41,1,5,12,15,16,29,36]; % well to process

% set experimental setting
frame_interval = 1; % hr
scale = 1.266; % um/px

%% for each well, estimate speed
for i = 1:length(pos)

    % load setting
    load("..\data\setting_"+pos(i)+".mat")
    
    % if there is a previous measurement, don't save
    if exist("..\data\kymograph_"+pos(i)+".mat",'file')
        isSave = 0;
    else
        isSave = 1;
    end

    % load images
    list = string(ls("..\img\*0.tif"));
    list = list(contains(list,pos(i)+"_") & contains(list,["cy5" "ros"])); % 2 signals are allowed
    img = cell(size(list));
    for j = 1:length(list)
        info = imfinfo("..\img\"+list(j));
        img{j} = arrayfun(@(i)imread("..\img\"+list(j), i, 'Info', info),...
            1:length(info),'UniformOutput',0);
    end

    % set direction for speed measurement
    K = kron([3 4 5;2 9 6;1 8 7],ones(ceil(size(initiation.mask)./3)));
    K = K(1:size(initiation.mask,1),1:size(initiation.mask,2));
    tbl = tabulate(K(initiation.mask));
    [~,index] = max(tbl(:,2));
    if index<9
        rotation_angle = mod((0:5:90)'+(45.*(index-1)),360);
    else
        rotation_angle = (0:5:359)';
    end

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

    % kymograph, for each conditions, get the kymograph
    [C1,C2]=meshgrid(1:length(list),1:length(rotation_angle));
    kymograph = table(list(C1(:)),rotation_angle(C2(:)),C1(:),C2(:),...
        'VariableNames',["source_img","angle","index_img","index_angle"]);

    for j = 1:height(kymograph)
        imgRotated = imrotate3(cat(3,img{kymograph.index_img(j)}{:}),...
            kymograph.angle(j),[0 0 1]);
        imgCropped = imgRotated(1:center(kymograph.index_angle(j),2),...
            center(kymograph.index_angle(j),1)-100:center(kymograph.index_angle(j),1)+100,:);
        maxR = max(imgCropped,[],2);
        kymograph.maxR{j} = rot90(maxR(index_strt(kymograph.index_angle(j)):end,:),-1);
    end

    % speed measurement
    % th = [0.1 0.3];
    th = [0.08 0.09]; % well 41
    % th = [0.1 0.15]; % well 16
    for j = 10%1:height(kymograph)
        maxR = kymograph.maxR{j};

        % remove all-black images if necessary
        index_frame = (1:size(maxR,1))'; % record the remaining frame
        index_noImg = find(sum(maxR==0,2)==size(maxR,2));
        if ~isempty(index_noImg)
            maxR(index_noImg,:) = [];
            index_frame(index_noImg) = [];
        end

        % identify the propagation front
        X = imbinarize(imadjust(maxR),th(kymograph.index_img(j)));
        X = imbinarize(imadjust(maxR),0.08);
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

        % skip non-monotonic change in distance
        if any(diff(x_mode)<-100)
            continue
        end

        % fit a straight line based on the curated data
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

    % check kymograph
    figure('WindowState','maximized')
    tiledlayout('flow','TileSpacing','compact','Padding','compact')
    order = reshape(1:height(kymograph),length(rotation_angle),[])';
    for j = order(:)'
        if ~isempty(kymograph.front{j}) % only plot when speed is measured
            nexttile
            hold on
            plotKymo(kymograph.maxR{j},frame_interval,scale)
            plot(kymograph.front{j}.time,kymograph.front{j}.mode_dist,'mv')
            hold off

            sig = extractBetween(kymograph.source_img(j),"_","_");
            title(sprintf("%s,%d,%.4f",sig,kymograph.angle(j),kymograph.speed{j}{1}(1)/3600))
            xticklabels("")
            yticklabels("")
        end
    end

    % save result
    if isSave
        exportgraphics(gcf,"..\img\check_speed_estimation\"+pos(i)+".jpg","Resolution",300) 
        save("..\data\kymograph_"+pos(i)+".mat","kymograph")
    end

    close all
end

%%
% % pos = "s"+[41,1,5,12,15,16,29,36]; % well to process
% clear
% load('..\data\kymograph_s36.mat')
% new = kymograph;
% load('..\data\v1\kymograph_ros_s36.mat')
% old = kymograph;
% load('..\data\v1\setting_s36.mat')
%
% index = find(ismember(0:5:359,speed.rotation_angle));
% [new.speed{index}{1}(1)/3600,old{1,5}{1}(1)/3600]
