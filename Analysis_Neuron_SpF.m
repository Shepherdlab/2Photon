function []=Analysis_Neuron_SpF(graph,Left,Right,Binoc)

% program is designed to import Suite2P data. Select F.npy, and then select
% the .xml data of your Prairie two-photon scan and .csv file of the Prairie
% two-photon electrophysiology. If anaylzing multiple recordings in a
% cocatenated stack, enter number of frames for each recording in the order
% in which they are in the stack. 

%graph: 1==display each ROI's traces, 0==only calculate and save outputs
%Left==# left frames, Right==# right frames, Binoc==# binoc frames

%Assumptions:
%In the same folder as your F.npy file, this program will look for your
%iscell.npy and Fneu.npy files that contain your manual classifications
%from Suite2P and the trace of the neuropil surrounding each neuron. In the
%same file as your Prairie data, this program will look for stimevents.m,
%which contains the identity of your stimuli associated with each 5 volt
%square wave pulse. stimevents is read left to right top to bottom. With
%each row being a repeated block of your randomized stimuli. In our case,
%this is a 5xNx2 matrix with N=#orientations x #spatial frequencies 

%Saves output in folder with Prairie data, separates output for each scan
%(ie, left, right, binoc). 

%Kyle Jenks, 2019-09-09. Shepherd Lab, University of Utah. 

%% open and import .npy file 
[filename, pathname] = uigetfile('F.npy', 'Select your .npy F data file');
%cancel if user clicks cancel
if isequal(filename,0) || isequal(pathname,0)
    disp('action canceled')
return
end
   %set current directory to pathname
   cd(pathname);
   %set file to path string
   file = [pathname filename]; 
   %Open Fluorescence info and neuropil data
   ROIdata_all = readNPY(file)';
  
   %open neuropil data
   Fneu = readNPY('Fneu.npy')';
   
   %open iscell
   iscell=readNPY('iscell.npy');
   %Only analyze neuropil and ROIs that user selected as cells (1=chosen,
   %0=discarded). 
   ROIdata_all=ROIdata_all(:,logical(iscell(:,1)));
   Fneu=Fneu(:,logical(iscell(:,1)));
   %Correction factor set to 100%. Fully subtract neuropil.
   ROIdata_sub = ROIdata_all-Fneu;
   
   %default number of loops
   number_loops=1;
   %If there are two scans
   if nargin==3
       if isequal(Left+Right,length(ROIdata_all))==0
           return
       end
   LeftFrames=ROIdata_sub(1:Left,:);
   RightFrames=ROIdata_sub(Left+1:Left+Right,:);
   number_loops=2;
   end
   %If there are three scans
    if nargin==4
       if isequal(Left+Right+Binoc,length(ROIdata_all))==0
           return
       end
   LeftFrames=ROIdata_sub(1:Left,:);
   RightFrames=ROIdata_sub(Left+1:Left+Right,:);
   BinocFrames=ROIdata_sub(Left+Right+1:Left+Right+Binoc,:);
   number_loops=3;
   end
   
   %% loop through left frames, then right, then binoc
   for hh=1:number_loops
       if hh==1
           if hh==1 && nargin>1
           ROIdata_raw=LeftFrames;
           else
               ROIdata_raw=ROIdata_sub;
           end
       end
       if hh==2
           ROIdata_raw=RightFrames;
       end
       if hh==3
           ROIdata_raw=BinocFrames;
       end
       
%% open and import XML file
  [filename, pathname] = uigetfile('*.xml', 'Select your .xml file from Prairie data');
if isequal(filename,0) || isequal(pathname,0)
    disp('action canceled')
    return
end
      %set current directory to pathname
    cd(pathname);
   file = [pathname filename]; 
   data = xml2struct(file);

%determine length
framelength = length(data.PVScan.Sequence.Frame);
%create vector of proper size for timestamps
timestamps=zeros(framelength,1);
%loop to extract timestamps, either relative or absolute
if (data.PVScan.Sequence.Frame{1, 2}.Attributes.relativeTime>0)
for i = 1:framelength
    timestamps(i,1)=str2double(data.PVScan.Sequence.Frame{1, i}.Attributes.relativeTime);
end
else
  for i = 1:framelength
    timestamps(i,1)=str2double(data.PVScan.Sequence.Frame{1, i}.Attributes.absoluteTime);
  end  
end
clear i

%% Select Electrophys data
%Electrophys file used to calculate timing of stimulation events
%open file
[filename, pathname] = uigetfile('*.csv', 'Select your .CSV file for electrophysiology');
if isequal(filename,0) || isequal(pathname,0)
    disp('action canceled')
    return
end
   file = [pathname filename]; 
   electrophystime = csvread(file,2);

%Shift trace up and down to find start and end of square wave pulses. 
roundedtimes=round(electrophystime);
shiftdown=[0,0;roundedtimes];
shiftup=[roundedtimes;0,0];
%find differences between the two shifted vectors
firstevent=(eq(shiftdown(:,2),shiftup(:,2)));
%1 now means there was an event
secondevent=1.-firstevent;
%event will now be same size as original roundedtimes matrix
thirdevent=secondevent(1:end-1);
%multiply together with rounded times column 1 to get timestamps where
%events happened. because it is relative time, this corresponds to
%handles.timestamps for image data (once divided by 1000 to be seconds).
event_timestamps=(roundedtimes(:,1).*thirdevent)/1000;
%remove zeros from event_timestamps 
Events_notzero=event_timestamps(event_timestamps~=0);
%remove events that follow the preceding event with less than 0.5 second
%delay. indicates an error introduced by a non-instantaneous change in
%voltage.
Events_notzero = Events_notzero(diff([0 Events_notzero'])>0.5);
%take every 2nd value of a matrix, this is the start of each block of
%stimuli (starting from 1st event)
startevent=Events_notzero(1:2:end);
%take every 2nd value from 2nd event, this is the blank frames start time
blankframes=Events_notzero(2:2:end);

%% find F for each ROI and convert data to Z-score
%calculate number of ROIs
[TimeLength,numROIs]=size(ROIdata_raw);
%preallocate matrix for blank frame data (period b/t stims, ISI)
F_data=zeros(1,(numROIs));
%loop to collect ISI (assumes FIVE seconds of inter stimulus interval)
for ii=1:length(blankframes)
    Chunk = (ROIdata_raw(timestamps(:,1) >=(blankframes(ii)+1)&timestamps(:,1) <=(blankframes(ii)+5),:));
    F_data=[F_data;Chunk];
end
%exclude first empty row
F_data=F_data(2:end,:);
clear ii 
%calculate mean and stdev of each ROI's ISI
meanF=mean(F_data(1:end,:),1);
stdevF=std(F_data(1:end,:),[],1);

%Z-Score data
ROIdata=zeros(TimeLength,numROIs);
for ii=1:numROIs
ROIdata(:,ii)=(ROIdata_raw(:,ii)-meanF(1,ii))./stdevF(1,ii);
end
clear ii 

%% process by event
%.mat file that has the stimulus order saved. 
load('stimevents');
[counts,~,third_dim] = size(stimevents);

if third_dim ==1
    print('Visual presentation was done using the older visual stimulation program. Revert to Analysis_Neuron')
    return
end

%find # of events
stim_orientations = stimevents(:,:,1)';
stim_frequencies = stimevents(:,:,2)';
stimnumber=numel(stim_orientations);

%Resphape into vectors and combine along with event times
stim_orientations = reshape(stim_orientations,1,stimnumber)';
stim_frequencies = reshape(stim_frequencies,1,stimnumber)';
stim_pairs = [stim_orientations,stim_frequencies,startevent];

%Sort by orientation
stim_pairs = sortrows(stim_pairs);

%Number of unique stims (orientation/frequency pairs) 
unique_orientations = sort(unique(stim_orientations));
unique_frequencies = sort(unique(stim_frequencies));

%% extract time series data. ASSUMES A FIVE SECOND STIM.
for ii = 1:length(unique_orientations)
    %find indexes of orientations
    orientation_current = find(stim_pairs(:,1)==unique_orientations(ii));
    for jj = 1:length(unique_frequencies)
        %find indexes of frequencies
        frequency_current = find(stim_pairs(:,2)==unique_frequencies(jj));
        %Find rows with the current orientation and frequency
        current_stim = intersect(orientation_current,frequency_current);
        for kk = 1:length(current_stim)
            Datacell{ii,jj,kk}=timestamps(timestamps >=(stim_pairs(current_stim(kk),3)-5)&timestamps <=(stim_pairs(current_stim(kk),3)+10));
       % ROI data for all ROIs for that time range around the stim
        Datacell{ii,jj,kk}(:,2:numROIs+1)=ROIdata(timestamps >=(stim_pairs(current_stim(kk),3)-5)&timestamps <=(stim_pairs(current_stim(kk),3)+10),:);
        %normalize times of each event to be on a scale with 0 as the event
        %by subtracting the stim time
        Datacell{ii,jj,kk}(:,1)=Datacell{ii,jj,kk}(:,1)-stim_pairs(current_stim(kk),3);   
        end
    end
end
clear ii jj kk orientation_current frequency_current 

%% Calculate mean response and p values
%Preallocate matrix for mean response to stimuli
mean_data=zeros(numROIs,length(unique_orientations),length(unique_frequencies));
%preallocate matrix for the t test of pre vs stim means
p_data=zeros(numROIs,length(unique_orientations),length(unique_frequencies));
%preallocate matrix for the mean pre-stim period
mean_pre=zeros(numROIs,length(unique_orientations),length(unique_frequencies));
%preallocate matrix for each individual pre trial period
pre_data=zeros(numROIs,length(unique_orientations),length(unique_frequencies),length(current_stim));
trial_data = zeros(numROIs,length(unique_orientations),length(unique_frequencies),length(current_stim));

%Calculate without graphing data
if nargin>0 && graph==0 || nargin==0
%Loop through the ROIs
for l=1:numROIs
    for ii = 1:length(unique_orientations)
    %find indexes of orientations
    orientation_current = find(stim_pairs(:,1)==unique_orientations(ii));
    for jj = 1:length(unique_frequencies)
        %find indexes of frequencies
        frequency_current = find(stim_pairs(:,2)==unique_frequencies(jj));
        %Find rows with the current orientation and frequency
        current_stim = intersect(orientation_current,frequency_current);
         %create empty matrix to contain interpolated stim traces
    stim_traces=zeros(151,length(current_stim));
        for kk = 1:length(current_stim)
       %interpolate each trial on a 0.1 S scale from -5 to 10 seconds 
       stim_traces(:,kk)=interp1(Datacell{ii,jj,kk}(:,1),Datacell{ii,jj,kk}(:,l+1),-5:0.1:10);
       %take mean of individual trials
       trial_data(l,ii,jj,kk)=mean(Datacell{ii,jj,kk}(Datacell{ii,jj,kk}(:,1)>=0&Datacell{ii,jj,kk}(:,1)<=5,l+1));
       %take mean of pre-stim periods from -4 to 0 seconds
       pre_data(l,ii,jj,kk)=mean(Datacell{ii,jj,kk}(Datacell{ii,jj,kk}(:,1)>=-4&Datacell{ii,jj,kk}(:,1)<=0,l+1));
        end
    %Take mean of traces using robust mean. This discards outliers for each
    %time point, which helps the mean curve to be more representative of an
    %actual response. Check robustMean documentation to determine if this
    %would be applicable to your analysis
    mean_trace=robustMean(stim_traces,2);
    %take the mean of the mean trace during stimulus period (0-5)
    mean_data(l,ii,jj)=mean(mean_trace(51:101,:));
    %take the mean of the mean trace during pre-stimulus period (-4 to 0)
    mean_pre(l,ii,jj)=mean(mean_trace(11:50,:));
    %statistical test to compare pre and stim means 
    [~,p_data(l,ii,jj)]=ttest(reshape(pre_data(l,ii,jj,:),length(current_stim),1),reshape(trial_data(l,ii,jj,:),length(current_stim),1));
clear mean_trace stim_traces TF 
    end
    end
end
end
clear ii jj kk

%Graph and calculate the data
if nargin>0 && graph==1
%Loop through the ROIs
for l=1:numROIs
    %create new figure with ROI name and position. Adjust for your monitor
    %needs. 
figure('name',sprintf('Plot of ROI %d',l),'numbertitle','off','position',[250,50,1000,700])
    %loop through stims
    for ii = 1:length(unique_orientations)
    %find indexes of orientations
    orientation_current = find(stim_pairs(:,1)==unique_orientations(ii));
    for jj = 1:length(unique_frequencies)
        %find indexes of frequencies
        frequency_current = find(stim_pairs(:,2)==unique_frequencies(jj));
        %Find rows with the current orientation and frequency
        current_stim = intersect(orientation_current,frequency_current);
         graphsize=length(unique_orientations)/2;     
       subplot(2,graphsize,ii);
    %create empty matrix to contain interpolated stim traces
    stim_traces=zeros(151,length(current_stim));
        for kk = 1:length(current_stim)
       %interpolate each trial on a 0.1 S scale from -5 to 10 seconds 
       stim_traces(:,kk)=interp1(Datacell{ii,jj,kk}(:,1),Datacell{ii,jj,kk}(:,l+1),-5:0.1:10);
       %take mean of individual trials
       trial_data(l,ii,jj,kk)=mean(Datacell{ii,jj,kk}(Datacell{ii,jj,kk}(:,1)>=0&Datacell{ii,jj,kk}(:,1)<=5,l+1));
       %take mean of pre-stim periods from -4 to 0 seconds
       pre_data(l,ii,jj,kk)=mean(Datacell{ii,jj,kk}(Datacell{ii,jj,kk}(:,1)>=-4&Datacell{ii,jj,kk}(:,1)<=0,l+1));
        end
    %Take mean of traces using robust mean. This discards outliers for each
    %time point, which helps the mean curve to be more representative of an
    %actual response. Check robustMean documentation to determine if this
    %would be applicable to your analysis
    mean_trace=robustMean(stim_traces,2);
    %plot mean trace
    plot(-5:.1:10,mean_trace,'LineWidth',2)
    hold on
    %take the mean of the mean trace during stimulus period (0-5)
    mean_data(l,ii,jj)=mean(mean_trace(51:101,:));
    %take the mean of the mean trace during pre-stimulus period (-4 to 0)
    mean_pre(l,ii,jj)=mean(mean_trace(11:50,:));
    %statistical test to compare pre and stim means 
     [~,p_data(l,ii,jj)]=ttest(reshape(pre_data(l,ii,jj,:),length(current_stim),1),reshape(trial_data(l,ii,jj,:),length(current_stim),1));
clear mean_trace stim_traces TF trial_data pre_data
    end
    legend(string((unique_frequencies)));
    hold off
    end
    %pause then advance. Allows user time to look at the plots. 
pause on
    pause(5);
pause off
close all
end
end
clear ii jj kk  

%Clear all unnecessary variables
 clear avgData avgH blankframes Chunk counts data databetween Datacell electrophystime event_timestamps eventmatrix Events_notzero F_data file 
 clear filename firstevent framelength graphsize i ii j l numROIs pathname ROIdata_raw roundedtimes secondevent shiftdown 
 clear shiftup standarddev startevent stim stimevents stimnumber thirdevent timestamps uniquestims uniquetimestamps w ROIdata ROIdata_raw ROIdata_all
 clear Fneu iscell ISI_FoverF Legend ROIdata_sub TimeLength
 %Save Data in the current directory
 save('Analysis_Neuron_SpF');
 %Report number of responsive cells in the command window
 Responsive_mean=mean_data>=0.5 & p_data<=0.05;
 Responsive_cells = max(Responsive_mean,[],3);
 Responsive_cells=sum(Responsive_cells,2);
 Responsive_cells(Responsive_cells>=1)=1;
 Responsive_cells=sum(Responsive_cells);
 Report=sprintf('There were %s responsive cells',num2str(Responsive_cells));
 disp(Report)
 clear Report Responsive_cells
 %Play a beep
 beep
 end 
end
