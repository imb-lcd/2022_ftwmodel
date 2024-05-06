%% initialize environment
clc
clear
close all

% available drug: DFO, FAC, Dasatinib, LY294002, GKT
% 24/96 wells
drug = "Dasatinib";
plate = 96;
% set experimental setting
frame_interval = 1; % hr
scale = 1.266; % um/px

% find related image filename
list = string(ls("..\img\"+drug+"_"+plate+"_*.tif"));
expression = ['(?<drug>\w+)_(?<plate>\d+)_s(?<well>\d+)_c(?<ch>\d+)_ff1_10.tif'];
tokenNames = regexp(list,expression,'names');
pos = struct2table([tokenNames{:}]');
pos.plate = str2double(pos.plate);
pos.well = str2double(pos.well);
pos.ch = str2double(pos.ch);

% alpha = [225 45 50 15];
% Dasatinib_24-2: 1:4770
% GKT_24-1: 1:5925
% GKT_24-2: 1:5090
% GKT_24-4: 1:5530

%% for each well, estimate speed
for i = 11%1:height(pos)    
    
    % load setting
    oname = "..\data\"+pos.drug(i)+"_"+pos.plate(i)+"\setting_s"+pos.well(i)+".mat";
    load(oname)
    
    % if there is a previous measurement, don't save
%     if exist("..\data\"+pos.drug(i)+"_"+pos.plate(i)+"\kymograph_s"+pos.well(i)+".mat",'file')
%         isSave = 1;
%     else        
%         isSave = 1;
%     end
    isSave = 0;
    
    % load images
    list = string(ls("..\img\"+pos.drug(i)+"_"+pos.plate(i)+"_s"+pos.well(i)+"_*.tif"));
    img = cell(size(list));
    for j = 1:length(list)
        info = imfinfo("..\img\"+list(j));
        img{j} = arrayfun(@(i)imread("..\img\"+list(j), i, 'Info', info),...
            1:length(info),'UniformOutput',0);
    end
    
    % set direction for speed measurement
    K = kron([3 4 5;2 9 6;1 8 7],ones(ceil(size(initiation.mask)./3)));
    tbl = tabulate(K(initiation.mask));
    [~,index] = max(tbl(:,2));
    if index<9
        rotation_angle = mod((0:5:90)+(45.*(index-1)),360);
    else
        rotation_angle = 0:5:359;
    end
    rotation_angle = 85;
    
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
    kymograph = table(list(C1(:)),rotation_angle(C2(:))',C1(:),C2(:),...
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
    for j = 1:height(kymograph)
        maxR = kymograph.maxR{j};%(:,1:5925);

        % remove all-black images if necessary
        index_frame = (1:size(maxR,1))'; % record the remaining frame
        index_noImg = find(sum(maxR==0,2)==size(maxR,2));
        if ~isempty(index_noImg)
            maxR(index_noImg,:) = [];
            index_frame(index_noImg) = [];
        end

        % identify the propagation front
        X = imbinarize(imadjust(maxR,stretchlim(maxR)),0.1);
        Y = X .* (1:size(X, 1)).';  % Replace 1s by the column index
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
                
        % if there is a change in speed
        if pos.plate(i)==24
            tbl = readtable("..\data\chemical_perturbation_image_summary.xlsx");
            tbl = tbl(tbl.plate==24&tbl.treatment==pos.drug(i)&tbl.well==pos.well(i),:);

            % if the information of treatment is no available, use ischange
            % to determine the frame 1 after treatment
            if ~isnan(tbl.frame1_treat)
                frame_strt = [frame_strt;tbl.frame1_treat+1];
                frame_stop = [tbl.frame1_treat-1;frame_stop];
            else
                TF = ischange(normalize(test,'range').*100,'linear','MaxNumChanges',3);
                frame_change = find(TF);
                frame_change(frame_change<frame_strt | frame_change>frame_stop) = [];
                if ~isempty(frame_change)
                    [~,index] = max(min([abs(frame_change-frame_strt) abs(frame_change-frame_stop)]));
                    frame_strt = [frame_strt;frame_change(index)+1];
                    frame_stop = [frame_change(index)-1;frame_stop];
                else
                    frame_strt = [frame_strt;frame_strt];
                    frame_stop = [frame_stop;frame_stop];
                end
            end
        end
% frame_stop(2) = 26
        
        front = cell(size(frame_strt,1),1);
        p = cell(size(frame_strt,1),1);

        for k = 1:size(frame_strt,1)
            front{k} = G(ismember(G.frame,frame_strt(k):frame_stop(k)),:);
            front{k}(cellfun('length',front{k}.px)==0,:) = [];

            % convert frames and pixels to real time and distance
            front{k}.time = (front{k}.frame-1).*frame_interval;
            front{k}.dist = cellfun(@(x)(x-1).*scale,front{k}.px,'UniformOutput',0);

            % find the neighbors around the mode of distance for each frame
            nNbr = 11;
            [f,xi] = cellfun(@(a)ksdensity(a),front{k}.dist,'UniformOutput',0);
            [~,I] = cellfun(@(a)max(a),f);
            x_mode = arrayfun(@(i)xi{i}(I(i)),1:length(I))';
            [~,I] = arrayfun(@(i)sort(abs(front{k}.dist{i}-x_mode(i))),1:length(I),'UniformOutput',0);
            front{k}.data_dist = arrayfun(@(i)front{k}.dist{i}(I{i}(1:min([nNbr length(front{k}.dist{i})]))),1:length(I),'UniformOutput',0)';
            front{k}.mode_dist = x_mode;

            % skip non-monotonic change in distance
            if any(diff(x_mode)<-200)
                continue
            end

            % fit a straight line based on the curated data
            data_t = repelem(front{k}.time,cellfun('length',front{k}.data_dist));
            data_x = [front{k}.data_dist{:}]';
            p{k} = polyfit(data_t,data_x,1);
        end

        % organize data
        kymograph.front{j} = front; 
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
            for k = 1:size(kymograph.frame_range{j},1)
                tmp = kymograph.front{j}{k};                
                plot(tmp.time,tmp.mode_dist,'v')

                if ~isempty(kymograph.speed{j}{k})
                    speed = kymograph.speed{j}{k}(1)/3600;
                    text(mean(tmp.time),mean(tmp.mode_dist),sprintf("%.4f",speed),...
                    'color','w','FontWeight','bold','HorizontalAlignment','right')
                end
            end
            hold off

            sig = "cy5";
            title(sprintf("%s,%d",sig,kymograph.angle(j)))
            xticklabels("")
            yticklabels("")
        end
    end
    
    % save result
    if isSave
        exportgraphics(gcf,"..\img\speed_check\estimation_"+...
            pos.drug(i)+"_"+pos.plate(i)+"_s"+pos.well(i)+".jpg",...
            "Resolution",300)        
        save(strrep(oname,"setting","kymograph"),"kymograph")
    end

    close all
end

%%
% pos = "s"+[19,4,8,9,10,21,22,23]; % well to process
% clear
% load('..\data\kymograph_s23.mat')
% new = kymograph;
% load('..\data\v1\kymograph_cy5_s23.mat')
% old = kymograph;
% load('..\data\v1\setting_s23.mat')
% 
% index = find(ismember(0:5:359,speed.rotation_angle));
% [new.speed{index}{1}(1)/3600,old{1,5}{1}(1)/3600]

% tmp = [];
% for i = 1:height(kymograph)
%     tmp = [tmp;min(diff(kymograph.front{i}.mode_dist))];
% end
% max(tmp)
