function [analysiscell,numgroups,numsubjects,groupnames,numtrials]=Group_Analysis(numgroups,numsubjects,numtrials)
%Code to analyze data from Analysis_Neuron or Analysis_Neuropil. 

%inputs
%numgroups = number of groups. 
%numsubjects = how many subjects per group. This should be a vector of
%individuals per group. For example, two groups of 6 would be [6,6].
%numtrials = number of trials per subject. Only two is supported right now,
%and assumes you are entering contra and ipsi data for each subject. 

%outputs
%analysiscell = stores all data from analysis for every animal, columns
%represent group and rows individual animals within each group.
%numgroups = number of groups. 
%numsubjects = how many subjects per group. This should be a vector of
%individuals per group. For example, two groups of 6 would be [6,6].
%groupnames = names of groups
%numtrials = number of trials per subject. Only two is supported right now,
%and assumes you are entering contra and ipsi data for each subject. 

%Kyle Jenks, 2019-04-18. Shepherd Lab, University of Utah. 

%% Iput information about the experiment
% How many groups are there?
if nargin<3
numgroups=inputdlg('How many groups?');
numgroups=str2double(numgroups{1});
end
 
%Number of individuals per group
if nargin<3
numsubjects = inputdlg( cellstr( num2str((1:numgroups).', 'Number of subjects group %d') ) );
numsubjects=str2double(numsubjects);
end
%Name of groups
groupnames = inputdlg( cellstr( num2str((1:numgroups).', 'Name of group %d') ) );

%Number of trials (currently only 2 supported)
if nargin<3
numtrials=inputdlg('How many trials per subject?');
numtrials=str2double(numtrials);
end

if numtrials==2  
%% Select all data for experiment
datacontra=cell(max(numsubjects),numgroups);
dataipsi=cell(max(numsubjects),numgroups);
for ii=1:numgroups
    %instructions for how to enter data for group 1
   waitfor(msgbox(num2str((ii).','select data for group %d. Select data for each subject in the order contra and then ipsi')));
    for jj=1:numsubjects(ii)

[filename, pathname] = uigetfile('*.mat', 'select contra data');
%cancel if user clicks cancel
if isequal(filename,0) || isequal(pathname,0)
    disp('action canceled by user')
    break
else
    %set current directory to pathname
    cd(pathname);
    %set file to path string
   file = [pathname filename]; 
   %Open Fluorescence info 
   datacontra{jj,ii} = load(file);
end  
[filename, pathname] = uigetfile('*.mat', 'Select ipsi data');
%cancel if user clicks cancel
if isequal(filename,0) || isequal(pathname,0)
    disp('action canceled by user')
    break
else
    %set current directory to pathname
    cd(pathname);
    %set file to path string
   file = [pathname filename]; 
   %Open Fluorescence info 
   dataipsi{jj,ii} = load(file);
end  
   
    end
end
    
%% Analysis of data
%create cell arrays to store data
analysiscell=cell(max(numsubjects),numgroups);

%loop through subjects and analyze data. Store each individuals data in
%ints own cell segment in analysiscell
for ii=1:numgroups
    for jj=1:numsubjects(ii)
        clear respcontra respipsi onlycontra only ipsi binoc contramax ipsimax sizematrix contra180 ipsi180
        %set negative numbers to 0
        datacontra{jj,ii}.mean_data(datacontra{jj,ii}.mean_data<0)=0;
         dataipsi{jj,ii}.mean_data(dataipsi{jj,ii}.mean_data<0)=0;
         
         %Find ROIs that are visually responsive to the contra and/or ipsi
         %eye
         [respcontra,~]=find(datacontra{jj,ii}.mean_data>=.5 & datacontra{jj,ii}.p_data<=0.05);
         [respipsi,~]=find(dataipsi{jj,ii}.mean_data>=.5 & dataipsi{jj,ii}.p_data<=0.05);
         
         %Sort ROIs
         respcontra=sort(respcontra);
         respipsi=sort(respipsi);
         
         %Identifiy ROIs responsive to one or both eyes
         onlycontra=setdiff(respcontra,respipsi,'sorted');
         onlyipsi=setdiff(respipsi,respcontra,'sorted');
         binoc=intersect(respcontra,respipsi);
         %Store their numbers
         analysiscell{jj,ii}.contra=onlycontra;
         analysiscell{jj,ii}.ipsi=onlyipsi;
         analysiscell{jj,ii}.binoc=binoc;
         %Store number of responsive ROIs that are contra, binoc or ipsi
         analysiscell{jj,ii}.numcontra=length(onlycontra);
         analysiscell{jj,ii}.numipsi=length(onlyipsi);
         analysiscell{jj,ii}.numbinoc=length(binoc);
         [analysiscell{jj,ii}.totalcells,~]=size(datacontra{jj,ii}.mean_data);
         
         %Store mean data of responsive ROIs
         analysiscell{jj,ii}.mean_contra=datacontra{jj,ii}.mean_data(onlycontra,:);
         analysiscell{jj,ii}.mean_binoccontra=datacontra{jj,ii}.mean_data(binoc,:);
         analysiscell{jj,ii}.mean_ipsi=dataipsi{jj,ii}.mean_data(onlyipsi,:);
         analysiscell{jj,ii}.mean_binocipsi=dataipsi{jj,ii}.mean_data(binoc,:);
                  
         %Store max peak pair for each ROI, and index
         [analysiscell{jj,ii}.contraMax(:,1),analysiscell{jj,ii}.contraMax(:,2)]=max(datacontra{jj,ii}.mean_data(onlycontra,:),[],2);
         [analysiscell{jj,ii}.ipsiMax(:,1),analysiscell{jj,ii}.ipsiMax(:,2)]=max(dataipsi{jj,ii}.mean_data(onlyipsi,:),[],2);
         [analysiscell{jj,ii}.BinocCMax(:,1),analysiscell{jj,ii}.BinocCMax(:,2)]=max(datacontra{jj,ii}.mean_data(binoc,:),[],2);
         [analysiscell{jj,ii}.BinocIMax(:,1),analysiscell{jj,ii}.BinocIMax(:,2)]=max(dataipsi{jj,ii}.mean_data(binoc,:),[],2);
         
         
         %Binocular ODI calculation
         [contramax(:,1),contramax(:,2)]=max(datacontra{jj,ii}.mean_data(binoc,:),[],2);
         [ipsimax(:,1),ipsimax(:,2)]=max(dataipsi{jj,ii}.mean_data(binoc,:),[],2);
          analysiscell{jj,ii}.ODI_B=zeros(analysiscell{jj,ii}.numbinoc,1);
         for kk=1:analysiscell{jj,ii}.numbinoc
            analysiscell{jj,ii}.ODI_B(kk)=ipsimax(kk,1)/(ipsimax(kk,1)+contramax(kk,1));
         end
         clear kk
         
         %Calculate correlation between contra and ipsi input to binocular
         %neurons
          for ll=1:length(binoc)
                     analysiscell{jj,ii}.binoccellcorr(ll)=corr(datacontra{jj,ii}.mean_data(binoc(ll),:)',dataipsi{jj,ii}.mean_data(binoc(ll),:)')';
          end
          
         %Angles for OSI and DSI calculations
         OSI_Angles=[0,0.5236,1.0472,1.5708,2.0944,2.618];
         DSI_Angles=[0,0.523598775598299,1.04719755119660,1.57079632679490,2.09439510239320,2.61799387799149,3.14159265358979,3.66519142918809,4.18879020478639,4.71238898038469,5.23598775598299,5.75958653158129];
         %find size of data matrix
         sizematrix=size(datacontra{jj,ii}.mean_data);
         %preallocate zeros matrices for speed
         analysiscell{jj,ii}.contra180=zeros(sizematrix(1,1),6);
         analysiscell{jj,ii}.ipsi180=zeros(sizematrix(1,1),6);
         %average on 180 cycle
        for kk=1:6
        analysiscell{jj,ii}.contra180(:,kk)=(datacontra{jj,ii}.mean_data(:,kk)+datacontra{jj,ii}.mean_data(:,kk+6))./2;
        end
        clear kk
        for kk=1:6
        analysiscell{jj,ii}.ipsi180(:,kk)=(dataipsi{jj,ii}.mean_data(:,kk)+dataipsi{jj,ii}.mean_data(:,kk+6))./2;
        end
        clear kk
        
        %Store 180 data for C-B-I cell types
        analysiscell{jj,ii}.contraonly180=analysiscell{jj,ii}.contra180(onlycontra,:);
        analysiscell{jj,ii}.binoccontra180=analysiscell{jj,ii}.contra180(binoc,:);
        analysiscell{jj,ii}.ipsionly180=analysiscell{jj,ii}.ipsi180(onlyipsi,:);
        analysiscell{jj,ii}.binocipsi180=analysiscell{jj,ii}.ipsi180(binoc,:);
        
        %contra calculations
        %preallocate for speed
        analysiscell{jj,ii}.pref_C=zeros(analysiscell{jj,ii}.numcontra,1);
        analysiscell{jj,ii}.OSI_C=zeros(analysiscell{jj,ii}.numcontra,1);
        analysiscell{jj,ii}.DSI_C=zeros(analysiscell{jj,ii}.numcontra,1);
        for kk=1:analysiscell{jj,ii}.numcontra
            %preferred angle.loops through contraonly ROIs and runs
            %calculations on those rows in contra180
            analysiscell{jj,ii}.pref_C(kk)=sum(analysiscell{jj,ii}.contra180(onlycontra(kk),:).*exp(2*1i*OSI_Angles),2);
            %OSI, for contra ROIs. loops through contraonly ROIs and runs
            %calculations on those rows in contra180
            analysiscell{jj,ii}.OSI_C(kk)=abs(analysiscell{jj,ii}.pref_C(kk)./sum(analysiscell{jj,ii}.contra180(onlycontra(kk),:),2));
            %DSI. loops through contraonly ROIs and runs
            %calculations on those rows in the mean data.
            analysiscell{jj,ii}.DSI_C(kk)=abs(sum(datacontra{jj,ii}.mean_data(onlycontra(kk),:).*exp(1i*DSI_Angles),2)./sum(datacontra{jj,ii}.mean_data(onlycontra(kk),:),2));
        end
         clear kk
     %convert preference to degrees, angle converts on a pi to -pi scale, so adjust back. 
            analysiscell{jj,ii}.pref_C=degrees(wrapTo2Pi(angle(analysiscell{jj,ii}.pref_C))./2);
            %scale is 0-180 direction, convert to orientation (ie 0 degree
             %direction is actually a 90 degree bar moving in the 0 degree
             %direction)
            for kk=1:numel(analysiscell{jj,ii}.pref_C)
                if analysiscell{jj,ii}.pref_C(kk)<=90
                    analysiscell{jj,ii}.pref_C(kk)=analysiscell{jj,ii}.pref_C(kk)+90;
                else
                    if analysiscell{jj,ii}.pref_C(kk)>90
                    analysiscell{jj,ii}.pref_C(kk)=analysiscell{jj,ii}.pref_C(kk)-90;
                    end
                end
            end
            
            %Scale to a 0-90 degree scale. 
            for kk=1:numel(analysiscell{jj,ii}.pref_C)
                    if analysiscell{jj,ii}.pref_C(kk)>90
                    analysiscell{jj,ii}.pref_C(kk)=abs(analysiscell{jj,ii}.pref_C(kk)-180);
                    end
            end
            
        %ipsi calculations
        %preallocate for speed
        analysiscell{jj,ii}.pref_I=zeros(analysiscell{jj,ii}.numipsi,1);
        analysiscell{jj,ii}.OSI_I=zeros(analysiscell{jj,ii}.numipsi,1);
        analysiscell{jj,ii}.DSI_I=zeros(analysiscell{jj,ii}.numipsi,1);
        for kk=1:analysiscell{jj,ii}.numipsi
            %preferred angle.loops through ipsionly ROIss and runs
            %calculations on those rows in ipsi180
            analysiscell{jj,ii}.pref_I(kk)=sum(analysiscell{jj,ii}.ipsi180(onlyipsi(kk),:).*exp(2*1i*OSI_Angles),2);
            %OSI, for ipsi ROIs. loops through ipsionly ROIs and runs
            %calculations on those rows in ipsi180
            analysiscell{jj,ii}.OSI_I(kk)=abs(analysiscell{jj,ii}.pref_I(kk)./sum(analysiscell{jj,ii}.ipsi180(onlyipsi(kk),:),2));
            %DSI. loops through ipsionly ROIs and runs
            %calculations on those rows in the mean data.
            analysiscell{jj,ii}.DSI_I(kk)=abs(sum(dataipsi{jj,ii}.mean_data(onlyipsi(kk),:).*exp(1i*DSI_Angles),2)./sum(dataipsi{jj,ii}.mean_data(onlyipsi(kk),:),2));
        end
         clear kk
        
            %convert preference to degrees, angle converts on a pi to -pi scale, so adjust back. 
            analysiscell{jj,ii}.pref_I=degrees(wrapTo2Pi(angle(analysiscell{jj,ii}.pref_I))./2);
            %scale is 0-180 direction, convert to 90-90 orientation
            %scale is 0-180 direction, convert to orientation (ie 0 degree
             %direction is actually a 90 degree bar moving in the 0 degree
             %direction)
            for kk=1:numel(analysiscell{jj,ii}.pref_I)
                if analysiscell{jj,ii}.pref_I(kk)<=90
                    analysiscell{jj,ii}.pref_I(kk)=analysiscell{jj,ii}.pref_I(kk)+90;
                else
                    if analysiscell{jj,ii}.pref_I(kk)>90
                    analysiscell{jj,ii}.pref_I(kk)=analysiscell{jj,ii}.pref_I(kk)-90;
                    end
                end
            end
            %Scale to a 0-90 degree scale. 
            for kk=1:numel(analysiscell{jj,ii}.pref_I)
                    if analysiscell{jj,ii}.pref_I(kk)>90
                    analysiscell{jj,ii}.pref_I(kk)=abs(analysiscell{jj,ii}.pref_I(kk)-180);
                    end
            end
            
        %binoc calculations
        %preallocate for speed
        analysiscell{jj,ii}.pref_B=zeros(analysiscell{jj,ii}.numbinoc,2);
        analysiscell{jj,ii}.OSI_B=zeros(analysiscell{jj,ii}.numbinoc,2);
        analysiscell{jj,ii}.DSI_B=zeros(analysiscell{jj,ii}.numbinoc,2);
        for kk=1:analysiscell{jj,ii}.numbinoc
            %preferred angle.loops through ipsionly ROIs and runs
            %calculations on those rows in ipsi180
            analysiscell{jj,ii}.pref_B(kk,1)=sum(analysiscell{jj,ii}.contra180(binoc(kk),:).*exp(2*1i*OSI_Angles),2);
            analysiscell{jj,ii}.pref_B(kk,2)=sum(analysiscell{jj,ii}.ipsi180(binoc(kk),:).*exp(2*1i*OSI_Angles),2);
            %OSI, for ipsi ROIs. loops through ipsionly ROIs and runs
            %calculations on those rows in ipsi180
            analysiscell{jj,ii}.OSI_B(kk,1)=abs(analysiscell{jj,ii}.pref_B(kk,1)./sum(analysiscell{jj,ii}.contra180(binoc(kk),:),2));
            analysiscell{jj,ii}.OSI_B(kk,2)=abs(analysiscell{jj,ii}.pref_B(kk,2)./sum(analysiscell{jj,ii}.ipsi180(binoc(kk),:),2));
            %DSI. loops through ipsionly rois and runs
            %calculations on those rows in the mean data.
            analysiscell{jj,ii}.DSI_B(kk,1)=abs(sum(datacontra{jj,ii}.mean_data(binoc(kk),:).*exp(1i*DSI_Angles),2)./sum(datacontra{jj,ii}.mean_data(binoc(kk),:),2));
            analysiscell{jj,ii}.DSI_B(kk,2)=abs(sum(dataipsi{jj,ii}.mean_data(binoc(kk),:).*exp(1i*DSI_Angles),2)./sum(dataipsi{jj,ii}.mean_data(binoc(kk),:),2));
        end
         clear kk
        %convert preference to degrees, angle converts on a pi to -pi scale, so adjust back. 
            analysiscell{jj,ii}.pref_B=degrees(wrapTo2Pi(angle(analysiscell{jj,ii}.pref_B))./2);
            %scale is 0-180 direction, convert to orientation (ie 0 degree
             %direction is actually a 90 degree bar moving in the 0 degree
             %direction)
            for kk=1:numel(analysiscell{jj,ii}.pref_B)
                if analysiscell{jj,ii}.pref_B(kk)<=90
                    analysiscell{jj,ii}.pref_B(kk)=analysiscell{jj,ii}.pref_B(kk)+90;
                else
                    if analysiscell{jj,ii}.pref_B(kk)>90
                    analysiscell{jj,ii}.pref_B(kk)=analysiscell{jj,ii}.pref_B(kk)-90;
                    end
                end
            end
            %Scale to a 0-90 degree scale. 
            for kk=1:numel(analysiscell{jj,ii}.pref_B)
                    if analysiscell{jj,ii}.pref_B(kk)>90
                    analysiscell{jj,ii}.pref_B(kk)=abs(analysiscell{jj,ii}.pref_B(kk)-180);
                    end
            end
            
    end
end
clear jj ii

end
end