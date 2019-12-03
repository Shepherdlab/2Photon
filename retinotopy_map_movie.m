%Variables relevant 
acq_rate = 5; %frame rate of camera in Hz
rep_rate = 0.1735769; %repetition rate of stimulus in Hz 
%rep_rate = 0.1056609; %repetition rate of stimulus in Hz 
repetitions = 40;

%Selectfile
[filename, pathname] = uigetfile('*.tif', 'Select Tif for analysis');
%cancel if user clicks cancel
if isequal(filename,0) || isequal(pathname,0)
    disp('action canceled by user')
    return
end
tic
    %set current directory to pathname
    cd(pathname);
    %set file to path string
   FileName = [pathname filename]; 
%get info for the image stack
Info = imfinfo(FileName);
%creat a blank matrix
mat_length = Info(1).Height;
mat_width = Info(1).Width;
mat_depth = length(Info);
imagestack = zeros(mat_length, mat_width, mat_depth);

%Open each frame and save as a Z plane of the matrix imageStack
for k = 1:mat_depth
    currentImage = imread(FileName, k, 'Info', Info);
    imagestack(:,:,k) = currentImage;
end 
clear Info

%calculate the slow response and subtract, double rep rate for filter
filter_length = 1/rep_rate*acq_rate*2;
stack_mov = movmean(imagestack,filter_length,3);
imagestack = imagestack - stack_mov;
clear stack_mov

%create deltaF movie by subtracting mean of first 10 seconds 
imagestack = imagestack - mean(imagestack(:,:,1:acq_rate*10),3);
%take off first 10 seconds (blank period)
imagestack = imagestack(:,:,acq_rate*10:end);

%Make an average response
%make the first block
mean_movie = imagestack(:,:,1:round(1/rep_rate*acq_rate));
%make a vector of rep start times
start_times = 1:round(1/rep_rate*acq_rate):repetitions*round(1/rep_rate*acq_rate);
for ii = 2:repetitions-1
    % segment of movie
    frame_segment = imagestack(:,:,start_times(ii):start_times(ii+1));
    %find which is shorter
    length_movie = min([size(mean_movie,3),size(frame_segment,3)]);
    %average only over shortest length
    mean_movie = (mean_movie(:,:,1:length_movie) + frame_segment(:,:,1:length_movie));   
end
%take mean to get average
mean_movie = mean_movie ./(repetitions-1);

% FFT
stack_fft = fft(mean_movie,[],3);

%normalize. Just to play movie.
min_image = min(mean_movie,[],[1,2,3]);
mean_movie = mean_movie - min_image;
max_image = max(mean_movie,[],[1,2,3]);
mean_movie = mean_movie ./max_image;

implay(mean_movie)

%save movie
outputFileName = 'mean_movie.tif';
for K=1:length(mean_movie(1, 1, :))
   imwrite(mean_movie(:, :, K), outputFileName, 'WriteMode', 'append',  'Compression','none');
end

clear imagestack

% Take first harmonic frequency, -1 to center on middle of monitor.
fft_map = -1*stack_fft(:,:,2);

Pyy = abs(fft_map);

toc

%plot magnitude 
figure('Name', 'Magnitude')
imagesc(Pyy,prctile(Pyy(:),[0.5,99.5]))
axis equal

%plot phase
figure('Name', 'Phase')
imagesc(wrapToPi(angle(fft_map)), [-pi pi])
colormap(gca, jet)
colorbar
axis equal

%Threshold Phase by magnitude.
figure('Name', 'Phase Thresholded')
mag_thresh = max(power_add,[],'all');
fft_thresh = power_add>=mag_thresh*.60; 
thresh_phase = wrapToPi(angle(fft_map));
h = imagesc(thresh_phase, [-pi pi]);
colormap(gca, jet)
colorbar
set(h,'AlphaData',fft_thresh)
axis equal

clear fft_thresh acq_rate currentImage filename FileName filter_length 
clear frame_segment h ii k length_movie mag_thresh mat_depth mat_length 
clear mat_width max_image mean_movie min_image pathname Pyy rep_rate
clear repetitions stack_fft start_times thresh_phase K outputFileName
