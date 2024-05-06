function plotKymo(X,frame_interval,scale)
%plotKymo plot kymograph
%   Detailed explanation goes here

if isempty(frame_interval)
    frame_interval = 1;
end
if isempty(scale)
    scale = 1;
end

% [~,I] = min(X(:));
% X(I) = 0;

len_px = max(size(X));
B = imresize(X,[len_px len_px]);
H = fspecial('average',13);
MotionBlur = imfilter(B,H,'replicate');
% MotionBlur = imadjust(uint16(MotionBlur));

imagesc([0-0.5 size(X,1)-1+0.5].*frame_interval,[0 size(X,2)-1].*scale,...
    MotionBlur')

xlim([0-0.5 size(X,1)-1+0.5].*frame_interval)
ylim([0 size(X,2)-1].*scale)
axis xy
colormap gray
    
end