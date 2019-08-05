%%  Script is to plot contra binoc and ipsi cells 
%This relies on the fall.mat output from Suite2P and the analysis from
%analysis_neuron for a contra and ipsi scan

%select Fall.mat file
[filename, pathname] = uigetfile('*.mat', 'select Fall.mat file');
 %set current directory to pathname
    cd(pathname);
    %set file to path string
   file = [pathname filename]; 
   %Open Fluorescence info 
  load(file);

    
%update iscell with the user curated output from iscell
iscell=readNPY('iscell.npy');

% Open contra and ipsi data
[filename, pathname] = uigetfile('*.mat', 'select contra data');
    %set current directory to pathname
    cd(pathname);
    %set file to path string
   file = [pathname filename]; 
   %Open Fluorescence info 
   datacontra{1,1} = load(file);
 
[filename, pathname] = uigetfile('*.mat', 'Select ipsi data');
    %set current directory to pathname
    cd(pathname);
    %set file to path string
   file = [pathname filename]; 
   %Open Fluorescence info 
   dataipsi{1,1} = load(file);

%Calculate contra, binoc, ipsi
 %set negative numbers to 0
        datacontra{1,1}.mean_data(datacontra{1,1}.mean_data<0)=0;
         dataipsi{1,1}.mean_data(dataipsi{1,1}.mean_data<0)=0;
         
         %Find ROIs that are visually responsive to the contra and/or ipsi
         %eye
         [respcontra,~]=find(datacontra{1,1}.mean_data>=.5 & datacontra{1,1}.p_data<=0.05);
         [respipsi,~]=find(dataipsi{1,1}.mean_data>=.5 & dataipsi{1,1}.p_data<=0.05);
         
         %Sort ROIs
         respcontra=sort(respcontra);
         respipsi=sort(respipsi);
         
         %Identifiy ROIs responsive to one or both eyes
         onlycontra=setdiff(respcontra,respipsi,'sorted');
         onlyipsi=setdiff(respipsi,respcontra,'sorted');
         binoc=intersect(respcontra,respipsi);
         
%Open mean image
[filename, pathname] = uigetfile('*.tif', 'select mean image');
 
   %set current directory to pathname
   cd(pathname);
   %set file to path string
   file = [pathname filename]; 
   %Open 
   image = imread(file);
   
%measure size of image
[ypixels,xpixels]=size(image);

%make blank ROI fields for each type of neuron
cROI=zeros(ypixels,xpixels);
bROI=zeros(ypixels,xpixels);
iROI=zeros(ypixels,xpixels);

%Make background colors color coded
bluebg=zeros(ypixels,xpixels,3);
bluebg(:,:,3)=1;
greenbg=zeros(ypixels,xpixels,3);
greenbg(:,:,2)=1;
yellowbg=zeros(ypixels,xpixels,3);
yellowbg(:,:,1)=1;
yellowbg(:,:,2)=1;

%Remove discarded ROIs from iscell
stat=stat(1,logical(iscell(:,1)));

%collect neuron ROI positions and weights for all ROIs in each category

%contra
for jj=1:length(onlycontra)   
for ii=1:length(stat{1,onlycontra(jj)}.lam)
cROI(stat{1,onlycontra(jj)}.ypix(ii),stat{1,onlycontra(jj)}.xpix(ii))=stat{1,onlycontra(jj)}.lam(ii);
end
end
%binoc
for jj=1:length(binoc)   
for ii=1:length(stat{1,binoc(jj)}.lam)
bROI(stat{1,binoc(jj)}.ypix(ii),stat{1,binoc(jj)}.xpix(ii))=stat{1,binoc(jj)}.lam(ii);
end
end
%ipsi
for jj=1:length(onlyipsi)   
for ii=1:length(stat{1,onlyipsi(jj)}.lam)
iROI(stat{1,onlyipsi(jj)}.ypix(ii),stat{1,onlyipsi(jj)}.xpix(ii))=stat{1,onlyipsi(jj)}.lam(ii);
end
end

%plot all 3 types of neurons 
figure,imshow(0.5*imadjust(image))
hold all
h=imshow(bluebg);
set(h,'AlphaData',(cROI)./(.2*max(max(cROI))))
h=imshow(greenbg);
set(h,'AlphaData',(bROI)./(.2*max(max(bROI))))
h=imshow(yellowbg);
set(h,'AlphaData',(iROI)./(.2*max(max(iROI))))

