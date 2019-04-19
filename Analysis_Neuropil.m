function []=Analysis_Neuropil(graph,Left,Right,Binoc)

% program is designed to import Suite2P data. Select Fneu.npy, and then select
% the .xml of your Prairie two-photon scan and .csv file of the Prairie
% two-photon electrophysiology. If anaylzing multiple recordings in a
% cocatenated stack, enter number of frames for each recording in the order
% in which they are in the stack. 

%graph: 1==display each ROI's traces, 0==only calculate and save outputs
%Left==# left frames, Right==# right frames, Binoc==# binoc frames

%Assumptions:
%In the same folder as your F.npy file, this program will look for your
%iscell.npy and Fneu.npy files that contain your manual classifications
%from Suite2P and the trace of the neuropil surrounding each neuron. 
%Additionally, a .csv file labeled blood vessel.csv should contain a column
%vector of the mean fluorescence of your chosen blood vessel from every
%frame of your cocatenated stack. THis can be done manually and saved from
%ImageJ using the ROI manager. 
%in the same file as your prairie data, this program will look for 
%stimevents.m, which contains the identity of your stimuli associated with
%each 5 volt square wave pulse. stimevents is read left to right top to 
%bottom. 

%Saves output in folder with Prairie data, separate output for ech scan
%(ie, left, right, binoc). 

%Kyle Jenks, 2019-04-18. Shepherd Lab, University of Utah. 

%% open and import .npy file of left and right data
[filename, pathname] = uigetfile('Fneu.npy', 'Select your .npy Fneu data file');
%cancel if user clicks cancel
if isequal(filename,0) || isequal(pathname,0)
    disp('action canceled')
return
end
    %set current directory to pathname
    cd(pathname);
    %set file to path string
   file = [pathname filename]; 
   %Open Fluorescence info and neurophil data
  % ROIdata_raw = csvread(file);
   ROIdata_all = readNPY(file)';
  
   %open blood vessel data
   Fvessel = csvread('blood vessel.csv');
   
   %subtract fvessel from every column
   ROIdata_all=ROIdata_all-Fvessel;
   
   %open iscell
   iscell=readNPY('iscell.npy');
   %Only analyze neuropil around user selected cells (1=chosen,
   %0=discarded). 
   ROIdata_sub=ROIdata_all(:,logical(iscell(:,1)));
   
  
  %default number of loops
   number_loops=1;
   if nargin==3
       if isequal(Left+Right,length(ROIdata_all))==0
           return
       end
   LeftFrames=ROIdata_sub(1:Left,:);
   RightFrames=ROIdata_sub(Left+1:Left+Right,:);
   number_loops=2;
   end
    if nargin==4
       if isequal(Left+Right+Binoc,length(ROIdata_all))==0
           return
       end
   LeftFrames=ROIdata_sub(1:Left,:);
   RightFrames=ROIdata_sub(Left+1:Left+Right,:);
   BinocFrames=ROIdata_sub(Left+Right+1:Left+Right+Binoc,:);
   number_loops=3;
   end
   
   %% loop through left frames, then right
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
%create cell of proper size
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
%find differences
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
%take every 2nd value of a matrix, this is the start of each block of
%stimuli
startevent=Events_notzero(1:2:end);
%take every 2nd value form 2, this is the blank frames
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
%calculate mean and stdev of each ROI ISI
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

%find # of events
stimevents=stimevents';
stimnumber=numel(stimevents);
%reshape into vector.
stimevents=reshape(stimevents,1,stimnumber);

%organize into matrices
eventmatrix(:,1)=stimevents;
eventmatrix(:,2)=startevent(1:length(stimevents));

%number of unique stims
uniquestims=unique(eventmatrix(:,1));
%presentations per stim 
counts = histc(eventmatrix(:,1), uniquestims);

%% extract time series data. ASSUMES A FIVE SECOND STIM.
for i=1:length(uniquestims)
    %timestamps for uniquestim i
    stim=uniquestims(i);
    uniquetimestamps=eventmatrix(eventmatrix(:,1)==stim,2);
    %cell array to store data
    Datacell(i,1)={stim};
    %loop to get all events for each unique stim
    for j=1:counts(i)
       %data time for each stim, 5 seconds before to 10 seconds after stim
       %start time
       Datacell{i,j+1}=timestamps(timestamps >=(uniquetimestamps(j)-5)&timestamps <=(uniquetimestamps(j)+10));
       % ROI data for all ROIs for that time range around the stimulus
        Datacell{i,j+1}(:,2:numROIs+1)=ROIdata(timestamps >=(uniquetimestamps(j)-5)&timestamps <=(uniquetimestamps(j)+10),:);
        %normalize times of each event to be on a scale with 0 as the event
        %by subtracting the stim time
        Datacell{i,j+1}(:,1)=Datacell{i,j+1}(:,1)-uniquetimestamps(j);
    end
end
 clear i j
 
%% Calculate mean response and p values
%Preallocate matrix for mean response to stimuli
mean_data=zeros(numROIs,length(uniquestims));
%preallocate matrix for the T test of pre vs stim means
p_data=zeros(numROIs,length(uniquestims));
%preallocate matrix for each inidivudal pre trial period
pre_data=zeros(numROIs,length(uniquestims),max(counts));
%mean pre-stim period
mean_pre=zeros(numROIs,length(uniquestims));
%preallocate matrix for mean of each individual trial
trial_data=zeros(numROIs,length(uniquestims),max(counts));
%preallocate matrix for each inidivudal pre trial period
pre_data=zeros(numROIs,length(uniquestims),max(counts));

%Calculate without graphing data
if nargin>0 && graph==0 || nargin==0
%Loop through the ROIs
for l=1:numROIs
    %loop through stims
for i=1:length(uniquestims)
    %create empty matrix to contain interpolated stim traces
    stim_traces=zeros(151,counts(i));
    %loop through the number of stim repeats
    for j=1:counts(i)
       %interpolate each trial on a 0.1 S scale from -5 to 10 seconds 
       stim_traces(:,j)=interp1(Datacell{i,j+1}(:,1),Datacell{i,j+1}(:,l+1),-5:0.1:10);
       %take mean of individual trials
       trial_data(l,i,j)=mean(Datacell{i,j+1}(Datacell{i,j+1}(:,1)>=0&Datacell{i,j+1}(:,1)<=5,l+1));
       %take mean of pre-stim periods from -4 to 0 seconds
       pre_data(l,i,j)=mean(Datacell{i,j+1}(Datacell{i,j+1}(:,1)>=-4&Datacell{i,j+1}(:,1)<=0,l+1));
    end
    %Take mean of traces
    mean_trace=robustMean(stim_traces,2);
    %take the mean of the mean trace during stimulus period (0-5)
    mean_data(l,i)=mean(mean_trace(51:101,:));
    %take the mean of the mean trace during pre-stimulus period (-4 to 0)
    mean_pre(l,i)=mean(mean_trace(11:50,:));
    %statistical test to compare pre and stim means 
    [~,p_data(l,i)]=ttest(reshape(pre_data(l,i,:),counts(i),1),reshape(trial_data(l,i,:),counts(i),1));
clear mean_trace stim_traces TF
end

end
end

%Graph and calculate the data
if nargin>0 && graph==1
%Loop through the ROIs
for l=1:numROIs
    %create new figure with ROI name
figure('name',sprintf('Plot of ROI %d',l),'numbertitle','off','position',[250,500,1000,700])
    %loop through stims
for i=1:length(uniquestims)
    %create empty matrix to contain interpolated stim traces
    stim_traces=zeros(151,counts(i));
    %loop through the number of stim repeats
    for j=1:counts(i)
       %interpolate each trial on a 0.1 S scale from -5 to 10 seconds 
       stim_traces(:,j)=interp1(Datacell{i,j+1}(:,1),Datacell{i,j+1}(:,l+1),-5:0.1:10);
       %take mean of individual trials (0-5 seconds)
       trial_data(l,i,j)=mean(Datacell{i,j+1}(Datacell{i,j+1}(:,1)>=0&Datacell{i,j+1}(:,1)<=5,l+1));
       %take mean of pre-stim periods from -4 to 0 seconds
       pre_data(l,i,j)=mean(Datacell{i,j+1}(Datacell{i,j+1}(:,1)>=-4&Datacell{i,j+1}(:,1)<=0,l+1));
       %graph each trace
       graphsize=length(uniquestims)/2;     
       subplot(2,graphsize,i);
       plot(-5:.1:10,stim_traces(:,j))
       hold on
    end
    %Take mean of traces
    mean_trace=robustMean(stim_traces,2);
    %plot mean trace
    plot(-5:.1:10,mean_trace,'k','LineWidth',2)
    %take the mean of the mean trace during stimulus period (0-5)
    mean_data(l,i)=mean(mean_trace(51:101,:));
    %statistical test to compare pre and stim means 
    [~,p_data(l,i)]=ttest(reshape(pre_data(l,i,:),counts(i),1),reshape(trial_data(l,i,:),counts(i),1));
    %Insert a legend
    Legend=legend(sprintf('T-Test=%0.3f',p_data(l,i)),sprintf('Mean=%0.3f',mean_data(l,i))); 
    %Color plots with STD larger than 0.5 && P<0.05
if mean_data(l,i)>=0.5 && p_data(l,i)<=0.05
    Legend.FontWeight='bold';
    set(gca,'Color',[0.4 0.8 0.8]);
end
    clear mean_trace stim_traces TF
end
%pause then advance
pause on
    pause(5);
pause off
close all
end

end
clear l i j
%Clear all unnecessary variables
    clear avgData avgH blankframes Chunk counts data databetween Datacell electrophystime event_timestamps eventmatrix Events_notzero F_data file 
    clear filename firstevent framelength graphsize i ii j l numROIs pathname ROIdata_raw roundedtimes secondevent shiftdown 
 clear shiftup standarddev startevent stim stimevents stimnumber thirdevent timestamps uniquestims uniquetimestamps w ROIdata ROIdata_raw ROIdata_all
 clear Fneu iscell ISI_FoverF Legend ROIdata_sub TimeLength
 %Save Data in the CD
 save('Analysis_Neuropil');
 %Report number of responsive neuropil patches 
 Responsive_mean=mean_data>=0.5 & p_data<=0.05;
 Responsive_patches=(Responsive_mean);
 Responsive_patches=sum(Responsive_patches,2);
 Responsive_patches(Responsive_patches>=1)=1;
 Responsive_patches=sum(Responsive_patches);
 Report=sprintf('There were %s responsive patches',num2str(Responsive_patches));
 disp(Report)
 clear Report Responsive_patches
 %Play a beep
 beep
   end
end
