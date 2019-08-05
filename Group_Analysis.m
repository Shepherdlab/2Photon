function [analysiscell,numgroups,numsubjects,groupnames,numtrials]=Group_Analysis(numgroups,numsubjects,numtrials)
%Code to analyze data from Analysis_Neuron or Analysis_Neuropil. 

%inputs
%numgroups = number of groups. 
%numsubjects = how many subjects per group. This should be a vector of
%individuals per group. For example, two groups of 6 would be [6,6].
%numtrials = number of trials per subject. Only two is supported right now,
%and assumes you are entering contra and ipsi data for each subject. 

%If inputs are not specified, the function will prompt the user for the
%information. 

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
if nargin<1
numgroups=inputdlg('How many groups?');
numgroups=str2double(numgroups{1});
end
 
%Number of individuals per group
if nargin<2
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
   %Open data
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
   %Open data
   dataipsi{jj,ii} = load(file);
end  
   
    end
end
    
%% Analysis of data
%create cell arrays to store data, subjects x group
analysiscell=cell(max(numsubjects),numgroups);

%loop through subjects and analyze data. Store each individual's data in
%its own place in analysiscell
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
                  
         %Store max peak for each ROI and the index
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
          
         %Angles for OSI and DSI calculations in the order they are stored
         %from analysis_neuron and analysis_neuropil. Keep in mind a unique
         %stim has a direction and orientation that are distinct 
         OSI_Angles=[1.5708,1.0472,0.5236,0,2.6180,2.0944];
         DSI_Angles=[3.1415,2.6179,2.0943,1.5707,1.0471,0.5236,0,5.7595,5.2359,4.7123,4.1887,3.6651];
         %find size of data matrix
         sizematrix=size(datacontra{jj,ii}.mean_data);
         %preallocate zeros matrices 
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
        
        %Store 180 data for contraonly, binoc, and ipsionly cell types.
        %Binoc ROIs have their contra and ipsi responses analyzed
        %separately
        analysiscell{jj,ii}.contraonly180=analysiscell{jj,ii}.contra180(onlycontra,:);
        analysiscell{jj,ii}.binoccontra180=analysiscell{jj,ii}.contra180(binoc,:);
        analysiscell{jj,ii}.ipsionly180=analysiscell{jj,ii}.ipsi180(onlyipsi,:);
        analysiscell{jj,ii}.binocipsi180=analysiscell{jj,ii}.ipsi180(binoc,:);
        
        %contra calculations
        %preallocate for speed
        analysiscell{jj,ii}.pref_C=zeros(analysiscell{jj,ii}.numcontra,1);
        analysiscell{jj,ii}.OSI_C=zeros(analysiscell{jj,ii}.numcontra,1);
        analysiscell{jj,ii}.DSI_C=zeros(analysiscell{jj,ii}.numcontra,1);
        analysiscell{jj,ii}.prefD_C=zeros(analysiscell{jj,ii}.numcontra,1);
        for kk=1:analysiscell{jj,ii}.numcontra
            %preferred angle .loops through contraonly ROIs and runs
            %calculations on those rows in contra180
            analysiscell{jj,ii}.pref_C(kk)=sum(analysiscell{jj,ii}.contra180(onlycontra(kk),:).*exp(2*1i*OSI_Angles),2);
            %OSI, for contra ROIs. loops through contraonly ROIs and runs
            %calculations on those rows in contra180
            analysiscell{jj,ii}.OSI_C(kk)=abs(analysiscell{jj,ii}.pref_C(kk)./sum(analysiscell{jj,ii}.contra180(onlycontra(kk),:),2));
            %Calculate preferred direction
            analysiscell{jj,ii}.prefD_C(kk)=sum(datacontra{jj,ii}.mean_data(onlycontra(kk),:).*exp(1i*DSI_Angles),2);
              %DSI. loops through contraonly ROIs and runs
            %calculations on those rows in the mean data
            analysiscell{jj,ii}.DSI_C(kk)=abs(analysiscell{jj,ii}.prefD_C(kk)./sum(datacontra{jj,ii}.mean_data(onlycontra(kk),:),2));
        end
         clear kk
     %convert preferences to degrees, angle converts on a pi to -pi scale, so adjust back. 
            analysiscell{jj,ii}.pref_C=rad2deg(wrapTo2Pi(angle(analysiscell{jj,ii}.pref_C))./2);
            analysiscell{jj,ii}.prefD_C=rad2deg(wrapTo2Pi(angle(analysiscell{jj,ii}.prefD_C)));
           
            
        %ipsi calculations
        %preallocate for speed
        analysiscell{jj,ii}.pref_I=zeros(analysiscell{jj,ii}.numipsi,1);
        analysiscell{jj,ii}.OSI_I=zeros(analysiscell{jj,ii}.numipsi,1);
        analysiscell{jj,ii}.DSI_I=zeros(analysiscell{jj,ii}.numipsi,1);
        analysiscell{jj,ii}.prefD_I=zeros(analysiscell{jj,ii}.numipsi,1);
        for kk=1:analysiscell{jj,ii}.numipsi
            %preferred angle. loops through ipsionly ROIs and runs
            %calculations on those rows in ipsi180
            analysiscell{jj,ii}.pref_I(kk)=sum(analysiscell{jj,ii}.ipsi180(onlyipsi(kk),:).*exp(2*1i*OSI_Angles),2);
            %OSI, for ipsi ROIs. loops through ipsionly ROIs and runs
            %calculations on those rows in ipsi180
            analysiscell{jj,ii}.OSI_I(kk)=abs(analysiscell{jj,ii}.pref_I(kk)./sum(analysiscell{jj,ii}.ipsi180(onlyipsi(kk),:),2));
            %Calculate preferred direction
            analysiscell{jj,ii}.prefD_I(kk)=sum(dataipsi{jj,ii}.mean_data(onlyipsi(kk),:).*exp(1i*DSI_Angles),2);
            %DSI. loops through contraonly ROIs and  runs
            %calculations on those rows in the mean  data.
            analysiscell{jj,ii}.DSI_I(kk)=abs(analysiscell{jj,ii}.prefD_I(kk)./sum(dataipsi{jj,ii}.mean_data(onlyipsi(kk),:),2));
        end
         clear kk
        
            %convert preference to degrees, angle converts on a pi to -pi scale, so adjust back. 
            analysiscell{jj,ii}.pref_I=rad2deg(wrapTo2Pi(angle(analysiscell{jj,ii}.pref_I))./2);
            analysiscell{jj,ii}.prefD_I=rad2deg(wrapTo2Pi(angle(analysiscell{jj,ii}.prefD_I)));
           
        %binoc calculations
        %preallocate for speed
        analysiscell{jj,ii}.pref_B=zeros(analysiscell{jj,ii}.numbinoc,2);
        analysiscell{jj,ii}.OSI_B=zeros(analysiscell{jj,ii}.numbinoc,2);
        analysiscell{jj,ii}.DSI_B=zeros(analysiscell{jj,ii}.numbinoc,2);
        analysiscell{jj,ii}.prefD_B=zeros(analysiscell{jj,ii}.numbinoc,2);
        for kk=1:analysiscell{jj,ii}.numbinoc
            %preferred angle. loops through ROIs and runs
            %calculations on those rows for contra and ipsi 180 data
            analysiscell{jj,ii}.pref_B(kk,1)=sum(analysiscell{jj,ii}.contra180(binoc(kk),:).*exp(2*1i*OSI_Angles),2);
            analysiscell{jj,ii}.pref_B(kk,2)=sum(analysiscell{jj,ii}.ipsi180(binoc(kk),:).*exp(2*1i*OSI_Angles),2);
            %OSI, for binoc ROIs. loops through  ROIs and runs
            %calculations on those rows
            analysiscell{jj,ii}.OSI_B(kk,1)=abs(analysiscell{jj,ii}.pref_B(kk,1)./sum(analysiscell{jj,ii}.contra180(binoc(kk),:),2));
            analysiscell{jj,ii}.OSI_B(kk,2)=abs(analysiscell{jj,ii}.pref_B(kk,2)./sum(analysiscell{jj,ii}.ipsi180(binoc(kk),:),2));
            %Calculate preferred direction
            analysiscell{jj,ii}.prefD_B(kk,1)=sum(datacontra{jj,ii}.mean_data(binoc(kk),:).*exp(1i*DSI_Angles),2);
            analysiscell{jj,ii}.prefD_B(kk,2)=sum(dataipsi{jj,ii}.mean_data(binoc(kk),:).*exp(1i*DSI_Angles),2);
            %DSI. loops through ROIs and  runs
            %calculations on those rows in the mean  data for both contra
            %and ipsi
            analysiscell{jj,ii}.DSI_B(kk,1)=abs(analysiscell{jj,ii}.prefD_B(kk,1)./sum(datacontra{jj,ii}.mean_data(binoc(kk),:),2));
            analysiscell{jj,ii}.DSI_B(kk,2)=abs(analysiscell{jj,ii}.prefD_B(kk,2)./sum(dataipsi{jj,ii}.mean_data(binoc(kk),:),2));
        end
         clear kk
        %convert preference to degrees, angle converts on a pi to -pi scale, so adjust back. 
            analysiscell{jj,ii}.pref_B=rad2deg(wrapTo2Pi(angle(analysiscell{jj,ii}.pref_B))./2);
            analysiscell{jj,ii}.prefD_B=rad2deg(wrapTo2Pi(angle(analysiscell{jj,ii}.prefD_B)));
        %Find binocular offset. First subtract the preference through the
        %two eyes.
            analysiscell{jj,ii}.offset_B = abs(analysiscell{jj,ii}.pref_B(:,1)- analysiscell{jj,ii}.pref_B(:,2));
        %We have numbers greater than 90 because the math doesn't realize 0 and 180 are equivalent. If a number is greater than 90, subtract from 180 to correct it.    
            for kk=1:length(analysiscell{jj,ii}.offset_B)
                if analysiscell{jj,ii}.offset_B(kk)>90
                    analysiscell{jj,ii}.offset_B(kk)=180-analysiscell{jj,ii}.offset_B(kk);
                end
            end  
           
    end
end
clear jj ii

end
end