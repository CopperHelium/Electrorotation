% convert all the avi files to mat and save with the name of the related
% tdms file.

clc 
close all;
clear all;

%read the info in the xls file
xlsname = 'Data.xlsx';
sheet = 1;
xlRange = 'E3:F44';
[avi_list,tdms_list] = xlsread(xlsname,sheet,xlRange);



%start from the current path
start_path=pwd;
%select the folder wanted
folder_name = uigetdir(start_path,'Select folder wanting to convert the avi:');
%set the path
cd(folder_name);

%run for all the avi files in this folder;
for i=1:length(tdms_list)
    %get tdms name
    tdms_name=char(tdms_list(i));
    tdms_name=[tdms_name(1:2) 'h' tdms_name(3:5) 'm' tdms_name(6:end) 's.mat'];
    
    %judge whether this tdms file exist
    if exist(tdms_name)~=0
    %get related video name and set its path
    avi_name=avi_list(i);
    avi_name=[num2str(avi_name) '.avi'];
    avi_path=strcat(folder_name,'\',avi_name);
    %read the objects
    obj = VideoReader(avi_path); 
    % get the height and width for the video
    vidHeight = obj.Height;
    vidWidth = obj.Width;
    %calculate the frames
    vidFrames = round(obj.Duration*obj.FrameRate);
    T = obj.Duration/vidFrames;   % Sampling period
    Fs = obj.FrameRate;  % Sampling frequency of video        
    t = (0:vidFrames-1)*T; % Time vector
    %% load each frame
    ss=zeros(vidHeight,vidWidth,1,vidFrames);
% read each frame and convert to gray data
for ii=1:vidFrames
    clc
    disp(['Loading video ' strcat(num2str(ii/vidFrames),'%')]);
    imageframe=readFrame(obj);
    ss(:,:,1,ii)=0.2989*imageframe(:,:,1)+ 0.5870 * imageframe(:,:,2) + 0.1140 * imageframe(:,:,3);
end

% reshape the gray image date from 4-D to 3-D
ss=reshape(ss,[vidHeight vidWidth vidFrames]);
file_out = strcat(tdms_name(1:end-4),'-VIDEO.mat');
save(file_out,'ss','vidFrames','vidHeight','vidWidth','T','t','Fs');
    end
end

