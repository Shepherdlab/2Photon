%% select P14
%Selectfile
[filename, pathname] = uigetfile('*.mat', 'Select P14 data');
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
   load(FileName);

%set map names
map_1 = map_90;
map_2 = map_270;
map_3 = map_0;
map_4 = map_180;
clear map_90 map_270 map_0 map_180

%combine magnitudes
power_add = abs(map_1) + abs(map_2) + abs(map_3) + abs(map_4);

%plot combined magnitudes
figure('Name','Magnitude')
imagesc(power_add,prctile(power_add(:),[0.50,99.5]))
axis equal
box off
set(gca,'xtick',[])
set(gca,'ytick',[])

%calculate vertical subtracted phase using complex conjugate
phase_sub = map_1 + conj(map_2);

%plot new phase map
figure('Name', 'vertical Phase')
imagesc(wrapToPi(angle(phase_sub)), [-pi pi])
colormap(gca, jet)
colorbar
axis equal
box off
set(gca,'xtick',[])
set(gca,'ytick',[])

%Threshold Phase by magnitude
mag_thresh = max(power_add,[],'all');
binoc_thresh = power_add>=mag_thresh*.60; 
figure('Name', 'Vertical Phase Thresholded')
thresh_phase = wrapToPi(angle(phase_sub));
h = imagesc(thresh_phase, [-pi pi]);
colormap(gca, jet)
colorbar
set(h,'AlphaData',binoc_thresh)
axis equal
box off
set(gca,'xtick',[])
set(gca,'ytick',[])

figure('Name', 'binoc_thresh')
binoc_thresh_P14 = double(binoc_thresh);
binoc_thresh_P14(binoc_thresh==1) = 0.6; 
binoc_thresh_P14(binoc_thresh==0) = 1; 
imshow(binoc_thresh_P14);


%clear everything except the binoc thresh map
clear filename FileName h mag_thresh map_1 map_2 map_3 map_4 pathname
clear phase_sub power_add thresh_phase

%% select P20
[filename, pathname] = uigetfile('*.mat', 'Select ipsi data');
%cancel if user clicks cancel
if isequal(filename,0) || isequal(pathname,0)
    disp('action canceled by user')
    return
end

    %set current directory to pathname
    cd(pathname);
    %set file to path string
   FileName = [pathname filename]; 
   load(FileName);

%set map names
map_1 = map_90;
map_2 = map_270;
map_3 = map_0;
map_4 = map_180;
clear map_90 map_270 map_0 map_180

%combine magnitudes
power_add = abs(map_1) + abs(map_2) + abs(map_3) + abs(map_4);

%plot combined magnitudes
figure('Name','P20 Magnitude')
imagesc(power_add,prctile(power_add(:),[0.50,99.5]))
axis equal
box off
set(gca,'xtick',[])
set(gca,'ytick',[])

%calculate vertical subtracted phase using complex conjugate
phase_sub = map_1 + conj(map_2);

%plot new phase map
figure('Name', 'Ipsi vertical Phase')
imagesc(wrapToPi(angle(phase_sub)), [-pi pi])
colormap(gca, jet)
colorbar
axis equal
box off
set(gca,'xtick',[])
set(gca,'ytick',[])

%Threshold Phase by magnitude
mag_thresh = max(power_add,[],'all');
P20_thresh = power_add>=mag_thresh*.60; 
figure('Name', 'Ipsi Vertical Phase Thresholded')
thresh_phase = wrapToPi(angle(phase_sub));
h = imagesc(thresh_phase, [-pi pi]);
colormap(gca, jet)
colorbar
set(h,'AlphaData',P20_thresh)
axis equal
box off
set(gca,'xtick',[])
set(gca,'ytick',[])

figure('Name', 'ipsi_thresh')
binoc_thresh_P20 = double(P20_thresh);
binoc_thresh_P20(P20_thresh==1) = 0.3; 
binoc_thresh_P20(P20_thresh==0) = 1; 
imshow(binoc_thresh_P20);


figure('Name', 'P14, P20 overlap')
overlap_thresh = binoc_thresh_P14 + binoc_thresh_P20;
overlap_thresh = overlap_thresh-1;
imshow(overlap_thresh);


clear all
toc