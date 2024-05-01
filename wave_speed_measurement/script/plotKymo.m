function plotKymo(X,frame_interval,scale)
%plotKymo a function to plot kymograph
%   The function plots the kymograph based on the real sampling interval
%   (h) and the real size (micrometer) without return anything.
%   X: kymograph
%   frame_interval: sampling interval
%   scale: pixel size (assume pixel width and height are the same)

% if there is no frame_interval provided, the default is 1 hr per frame
if isempty(frame_interval)
    frame_interval = 1;
end

% if there in no scale provided, the default is 1 um per pixel
if isempty(scale)
    scale = 1;
end

% interpolate the kymograph
len_px = max(size(X)); % length in pixel
B = imresize(X,[len_px len_px]);
H = fspecial('average',11);
MotionBlur = imfilter(B,H,'replicate'); % Gaussian blue
MotionBlur = imadjust(uint16(MotionBlur)); % convert the kymograph to 16-bit

% plot the kymograp based on sampling interval and pixel size
imagesc([0-0.5 size(X,1)-1+0.5].*frame_interval,[0 size(X,2)-1].*scale,...
    MotionBlur')
xlim([0-0.5 size(X,1)-1+0.5].*frame_interval)
ylim([0 size(X,2)-1].*scale)
axis xy
colormap gray

end