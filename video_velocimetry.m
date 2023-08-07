%%
clear
close all

% Define the video file, start frame, and end frame
vidObj = VideoReader('IMG_4429.MOV');
start_frame = 1;
end_frame = vidObj.NumFrames;

% Display the video first frame and let the user click two points to define
% bottom and top of conduit.
vidFrame = readFrame(vidObj);
figure();
imshow(vidFrame);
[x,y] = getline();
% determine the orientation of the conduit
dx = x(2)-x(1);
dy = y(2)-y(1);
theta = 90-atand(dy/dx);
close();
% display the rotated slice and have the user select the column (vertical
% slice position) to be used in the analysis.
figure();
tmp = imrotate(vidFrame,-theta,'bilinear','crop');
imshow(tmp);
[x,y] = ginput(1);
col = round(x);
close();
%%
video_slice = zeros(vidObj.Height,vidObj.NumFrames,3,'uint8');
video_slice(:,1,:) = tmp(:,col,:);
%
i=2;
figure();
while hasFrame(vidObj)
    vidFrame = readFrame(vidObj);                           % get the frame
    tmp = imrotate(vidFrame(:,col-100:col+100,:),-theta,'bilinear','crop');      % rotate to vertical
    video_slice(:,i,:) = tmp(:,101,:);                 % extract the slice
    if ~mod(i,1000)                                         % display output every so often
        clf;
        imshow(tmp)
        title([num2str(i) '/' num2str(vidObj.NumFrames)])
        drawnow();
    end
    i = i + 1;
end

%% Show the streak image
figure()
imshow(video_slice );
% select the region for velocimetry
[x,y] = ginput(2);
%%
start_ind = fix(x(1));
end_ind = fix(x(2));
imroi = video_slice(:,start_ind:end_ind,:);
imgray = rgb2gray(imroi);
figure
imagesc(imgray);
colorbar;

ntot = size(imroi,2);
nskip=1;
% ncorr=10; % number of rows from the video to use in cross-correlation analysis.
corr_length = 120;% number of frames for correlation in each direction
corr_rows = fix(linspace(1300,1800,10)); % rows to use in correlation

corr_range = (corr_length+1):nskip:(ntot-corr_length-1);
rs = zeros(corr_length*2+1,length(corr_range));
peak_lag = zeros(length(corr_range),1);

ind=1;
for i=corr_range
    for j=2:length(corr_rows)
        % determine number of rows for which to compute the cross-correlation
        % corr_rows = fix(linspace(1,vidObj.height,ncorr));
        corr_region = imgray(corr_rows,i-corr_length:i+corr_length);
        [r,lags] = xcorr(detrend(double(corr_region(1,:))),detrend(double(corr_region(2,:))),corr_length,'normalized');
        rs(:,ind) = r;
        [rmax,i1] = max(r);
        peak_lag(ind) = lags(i1);
        ind = ind+1;
        % compute the cross-correlation between each pair of rows for lags of
    end
end

figure;
subplot(2,1,1);
pcolor(1:size(rs,2),lags,rs); shading flat; colorbar
subplot(2,1,2);
plot(1:size(rs,2),peak_lag)
%% second approach, based on correlating adjacent frames
corr_length = 400;
rs = zeros(size(corr_range));
peak_lag = zeros(size(corr_range));
ind=1;
for i = corr_range
    % correlate this frame with the next one
    corr_region = imgray(corr_rows(1):corr_rows(end),i:i+1);
    [r,lags] = xcorr( detrend(double(corr_region(:,1))),detrend(double(corr_region(:,2))),corr_length,'normalized');
    [peaks,peak_ind] = findpeaks(r,'MinPeakProminence',0.05);
    peak_lags = lags(peak_ind);
    i1 = find(peak_lags < 0 ,1,'last'); % first peak to the left of 0 corresponding to upward propagating signal       
    if ~isempty(i1)
        peak_lag(ind) = peak_lags(i1);
    else
        peak_lag(ind) = NaN;
    end
    ind = ind + 1;
end



% video_bw = rgb2gray(video_slice);
% [gmag,angle] = imgradient(video_bw);
% figure, imshow(gmag);
%%
