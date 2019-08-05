function [AUC_data,AUC_data01,stdpeak_datathreetimes,stdpeak_data100]=Analysis_Spontaneous(eye_specific,Contra,Ipsi)
%% Analyze spontaneous events during the interstimulis interval
%The purpose of this function is to extract spontaneous calcium activity
%during the grey screeen period between stimuli, excluding the first second
%post stimulus

%Eye_specific 1 = analyze contra recording (1st scan in cocatenated stack)
%2 = analyze ipsi recording
%Contra = number of contra frames
%Ipsi = number of ipsi frames


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
   %Open Fluorescence info and neuropil data1
   ROIdata_all = readNPY(file)';
  
   %open neuropil data
   Fneu = readNPY('Fneu.npy')';
   
   %open iscell
   iscell=readNPY('iscell.npy');
   %Only analyze neuropil and ROIs that user selected as cells (1=chosen,
   %0=discarded). 
   ROIdata_all=ROIdata_all(:,logical(iscell(:,1)));
   Fneu=Fneu(:,logical(iscell(:,1)));
   %Correction factor set to 0.7
   ROIdata_sub = ROIdata_all-0.7*Fneu;
   
   %analyze which eye? assumes F and Fneu are stack of contra and ipsi
   
   %exit if contra+ipsi is not the length of the F data file
   if isequal(Contra+Ipsi,length(ROIdata_all))==0
           return
   end
   %if eye specific is 1, analyze contra
   if eye_specific==1
   ROIdata_sub=ROIdata_sub(1:Contra,:);
   end
   %if eye specific is 2, analyze ipsi
   if eye_specific==2
   ROIdata_sub=ROIdata_sub(Contra+1:Contra+Ipsi,:);
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
%take every 2nd value from 2, this is the blank frames
blankframes=Events_notzero(2:2:end);

%% find F for each ROI and convert data to delta F/F
%calculate number of ROIs
[TimeLength,numROIs]=size(ROIdata_sub);
%preallocate matrix for blank frame data (period b/t stims, ISI)
F_data=zeros(1,(numROIs));
%loop to collect ISI (assumes FIVE seconds of inter stimulus interval, take 1-5 seconds)
for ii=1:length(blankframes)
    Chunk = (ROIdata_sub(timestamps(:,1) >=(blankframes(ii)+1)&timestamps(:,1) <=(blankframes(ii)+5),:));
    F_data=[F_data;Chunk];
end
%exclude first empty row
F_data=F_data(2:end,:);
clear ii Chunk
%calculate mean of each ROI ISI
meanF=mean(F_data(1:end,:),1);

%Report ROIs with means <=0.
negative_means=meanF<=0;

%Delta F/F the whole trace
DeltaData=zeros(TimeLength,numROIs);
for ii=1:numROIs
DeltaData(:,ii)=(ROIdata_sub(:,ii)-meanF(1,ii))./meanF(1,ii);
end
clear ii    

%% find F for each ROI and convert data to delta F/F
%preallocate matrix for blank frame data (period b/t stims, ISI)
F_deltadata=zeros(1,(numROIs));
%loop to collect ISI (assumes FIVE seconds of inter stimulus interval, take 1-5 seconds)
for ii=1:length(blankframes)
    Chunk = (DeltaData(timestamps(:,1) >=(blankframes(ii)+1)&timestamps(:,1) <=(blankframes(ii)+5),:));
    F_deltadata=[F_deltadata;Chunk];
end
%exclude first empty row
F_deltadata=F_deltadata(2:end,:);
clear ii Chunk
%take stdev of blank frames of delta trace
STDEVdeltaF=std(F_deltadata);

%% Measure spontaneous events

% Positive AUC
%retake those 4 second periods from the now Delta F/F converted traces
%preallocate cell to store traces for every ROI
Delta_ISI=cell(length(blankframes),1);

%eliminate negative delta F/F values so we count only positive deflections
DeltaData(DeltaData<0)=0;

%loop to collect ISI (assumes FIVE seconds of inter stimulus interval, take 1-5 seconds)
for ii=1:length(blankframes)
    Chunk = (DeltaData(timestamps(:,1) >=(blankframes(ii)+1)&timestamps(:,1) <=(blankframes(ii)+5),:));
    Delta_ISI{ii,1}=Chunk;
end
clear Chunk ii

%take AUC for every 4 second period of every ROI, spacing is 0.378796 seconds
AUC_data=zeros(numROIs,length(blankframes));
for ii=1:numROIs
    for jj=1:length(blankframes)
        %AUC of data using trapezoid method
        AUC_data(ii,jj)=trapz(0.378796,Delta_ISI{jj}(:,ii));
    end
end
clear ii jj 
%Average by ROI
AUC_data=mean(AUC_data,2);

% AUC from Keck et al., Neuron 2013.
%retake those 4 second periods from the now Delta F/F converted traces
%preallocate cell to store traces for every ROI
Delta_ISI01=cell(length(blankframes),1);

%Keep only Delta values above 10% Delta F over F
DeltaData01=DeltaData-0.1;
%eliminate negative delta F/F values so we count only positive deflections
DeltaData01(DeltaData01<0)=0;

%loop to collect ISI (assumes FIVE seconds of inter stimulus interval, take 1-5 seconds)
for ii=1:length(blankframes)
    Chunk = (DeltaData01(timestamps(:,1) >=(blankframes(ii)+1)&timestamps(:,1) <=(blankframes(ii)+5),:));
    Delta_ISI01{ii,1}=Chunk;
end
clear Chunk ii

%take AUC for every 4 second period of every ROI, spacing is 0.378796 seconds
AUC_data01=zeros(numROIs,length(blankframes));
for ii=1:numROIs
    for jj=1:length(blankframes)
        %AUC of data using trapezoid method
        AUC_data01(ii,jj)=trapz(0.378796,Delta_ISI01{jj}(:,ii));
    end
end
clear ii jj
AUC_data01=mean(AUC_data01,2);

% Find peaks
%preallocate blank matrices
stdpeak_datathreetimes=zeros(numROIs,length(blankframes));
stdpeak_data100=zeros(numROIs,length(blankframes));
%Turn warning off for no peaks found
warning('off',['signal:findpeaks:largeMinPeakHeight'])
for ii=1:numROIs
    for jj=1:length(blankframes)
        %number of peaks above 3* the standard deviation 
        stdpeak_datathreetimes(ii,jj)=numel(findpeaks(Delta_ISI{jj}(:,ii),'MinPeakHeight',3*STDEVdeltaF(ii)));
        %number of peaks above 1 (100% standard deviation), Keck et al.
        %2013 metric
        stdpeak_data100(ii,jj)=numel(findpeaks(Delta_ISI{jj}(:,ii),'MinPeakHeight',1));
    end
end
%Turn warning back on so we don't mess with other programs
warning('on',['signal:findpeaks:largeMinPeakHeight'])
%Take average by ROIs
stdpeak_datathreetimes=mean(stdpeak_datathreetimes,2);
stdpeak_data100=mean(stdpeak_data100,2);
clear ii jj

save('spontaneous','AUC_data','AUC_data01','stdpeak_data100','stdpeak_datathreetimes');

