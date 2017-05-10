% Run wavelet for the avi files('*-VIDEO.mat')

clc 
close all;
clear all;

%start from the current path
start_path=pwd;
%select the '*-VIDEO.mat' wanted
Video_name = uigetfile(start_path,'Select the *-VIDEO.mat to make WT');
Video_path=strcat(start_path,'\',Video_name);
load(Video_path);


%% tracking the center for the video
pos=zeros(vidFrames,2);
for i=1:vidFrames
clc
disp(['Processing video ' strcat(num2str(i/vidFrames),'%')])
pic=ss(:,:,i);
pos(i,:)=center03(pic);
end

%% apply wavelet
wlet = 'morse';              % name of the mother wavelet
Scale_min=0;  % Set the minimum of frequency to 0 first, and search a proper one in the following steps
Scale_max=1e7; % Set the maximum of frequency to a large value first, and search a proper one in the following steps
    
% Got abs of coefficients for both X signals and Y signals. Plus them.
[coefsx, freqx] = cwt(pos(:,1),wlet,Fs,'TimeBandwidth',60,'VoicesPerOctave',48); 
[coefsy, freqy] = cwt(pos(:,2),wlet,Fs,'TimeBandwidth',60,'VoicesPerOctave',48);
Scale_min=max(Scale_min,find(freqx>20, 1, 'last' ));
Scale_max=min(Scale_max,length(freqx));
coepow=abs(coefsy)+abs(coefsx);
coepow(Scale_max+1:end,:)=[];
coepow(1:Scale_min-1,:)=[];
freq =freqx(Scale_min:1:Scale_max);

file_out = strcat(Video_name(1:end-10),'-AVIWAV.mat');
save(file_out,'pos','coepow','freq','t');
