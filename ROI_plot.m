function ROI_plot(ROI)

%%  Script is to plot ROI(s) of interest
%ROI should be scalar or vector of ROIs
%This relies on the fall.mat output from Suite2P. Mean image of imaged
%region should be a tif file with the same X by Y pixels as the data you
%analyzed. 

%select Fall.mat file
[filename, pathname] = uigetfile('*.mat', 'select Fall.mat file');
 %set current directory to pathname
    cd(pathname);
    %set file to path string
   file = [pathname filename]; 
   %Open Fluorescence info 
  load(file);

    
%update iscell in Fall.mat with the user curated output from iscell
iscell=readNPY('iscell.npy');
         
%Open image
[filename, pathname] = uigetfile('*.tif', 'select mean image');
 
    %set current directory to pathname
    cd(pathname);
    %set file to path string
   file = [pathname filename]; 
   %Open Fluorescence info 
   image = imread(file);

%measure size of image
[ypixels,xpixels]=size(image);

%make blank ROI fields for each type of neuron
cROI=zeros(ypixels,xpixels);

%Make background color blue
bluebg=zeros(ypixels,xpixels,3);
bluebg(:,:,3)=1;

%Remove discarded ROIs from stat based on iscell
stat=stat(1,logical(iscell(:,1)));

%collect ROI positions and pixel weights for all ROIs and place in your
%blue background
for jj=1:length(ROI)   
for ii=1:length(stat{1,ROI(jj)}.lam)
cROI(stat{1,ROI(jj)}.ypix(ii),stat{1,ROI(jj)}.xpix(ii))=stat{1,ROI(jj)}.lam(ii);
end
end

%Show the mean image, adjusted.
figure,imshow(0.5*imadjust(image))
hold all
%Now plot the blue background
h=imshow(bluebg);
%Set the transparency of the blue layer based on the values from the ROI
%weights, scaling to keep them in range and not labeling extraneuous pixels
set(h,'AlphaData',(cROI)./(.2*max(max(cROI))))



