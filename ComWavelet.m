% get the wavecoefficients for combined data

clc 
close all;
clear all;

% select '*s.mat'
start_path=pwd;
data_chose = uigetfile(start_path,'Select the *s.mat');
data_path=strcat(start_path,'\',data_chose);
load(data_path);
%load saved '-AVIMAV.mat' data
Video_name=[data_chose(1:end-4) '-AVIWAV.mat'];
load(Video_name);
coepow_avi=coepow;
coepow_avi=coepow_avi./(max(max(coepow_avi)));
freq_avi=freq;
avi_timevector=t; 
clear coepow;
clear freq;
clear t;

%load saved '-TDMSWAV.mat' data
TDMS_name=[data_chose(1:end-4) '-TDMSWAV.mat'];
load(TDMS_name);
coepow_tdms=coepow;
coepow_tdms=coepow_tdms./max(max(coepow_tdms));
freq_tdms=freq;
tdms_subsamp=sub_sampling;
tdms_timevector=t;
clear coepow;
clear freq;
clear t;

%% Set OFF and REC pairs for tdms files
    X_signal = UntitledVoltage_0.Data; 
    Y_signal = UntitledVoltage_1.Data; 
    fs = UntitledVoltage_2.Property.wf_samples; 
    F_state = UntitledVoltage_2.Data;
    F_state = (F_state > 2); 
    Switch = diff(F_state); 
    OnOff = find(Switch);  
    ON_tdms_pairs = [OnOff(1:2:end) OnOff(2:2:end)];
    OFF_tdms_pairs = [OnOff(2:2:end-2) OnOff(3:2:end)];
    numoff=size(OFF_tdms_pairs,1);
    if numoff>=3
        numoff=3;
    end
    numon=size(ON_tdms_pairs,1);
    if numon>=3
        numon=3;
    end
    ON_cnt = (ON_tdms_pairs(numon,2)-ON_tdms_pairs(numon,1)); % # samples in each ON event
    OFF_cnt = (OFF_tdms_pairs(numoff,2)-OFF_tdms_pairs(numoff,1)); % # samples in each OFF event 
    rec_flag = 0;
    if size(F_state,1)-ON_tdms_pairs(end,2)>OFF_cnt && F_state(end)==0
        OFF_tdms_pairs = vertcat([1 OnOff(1)],OFF_tdms_pairs,[OnOff(end) OnOff(end)+(OFF_tdms_pairs(numoff,2)-OFF_tdms_pairs(numoff,1))]);
        rec_flag = 1;
    elseif size(F_state,1)-ON_tdms_pairs(end,2)<OFF_cnt && F_state(end)==0
        OFF_tdms_pairs = vertcat([1 OnOff(1)],OFF_tdms_pairs);
        elsem
        OFF_tdms_pairs = vertcat([1 OnOff(1)],OFF_tdms_pairs);
    end
     REC_array = OFF_tdms_pairs(end,2):OFF_cnt:size(X_signal,1);
     REC_tdms_pairs = [REC_array(1:end-1)' REC_array(2:end)'];
     OFF_periods= (OFF_tdms_pairs-1)./fs;
     REC_periods=(REC_tdms_pairs-1)./fs;

     
     %% get info from '-VIDEO.mat' file
VID=[data_chose(1:end-4) '-VIDEO.mat'];
OUT=load(VID,'Fs','t');
Fs=OUT.Fs;
t=OUT.t;

%% Find OFF and REC pairs by related tdms files
% determine offset first
value_min=0.04; 
pre_peak=wltpeakfinder(coepow_avi(:,1:20*Fs),3,value_min);
pre_peak(isnan(pre_peak))=0;
pre_switch=diff(pre_peak);
off_start=find(abs(pre_switch)>20,1,'first');
if off_start<100 
    off_st=find(abs(pre_switch)>30,2,'first');
    off_start=off_st(2);
end
OffSet=off_start-floor(OFF_periods(1,2)*Fs)-1;
% set off pairs for avi
OFF_avi_pairs=zeros(size(OFF_tdms_pairs,1),2); % OFF pairs in the video 
for i=1: size(OFF_tdms_pairs,1)
OFF_avi_pairs(i,1) = find((t>=OFF_periods(i,1)),1,'first')+OffSet; % find the start point for OFF period i
OFF_avi_pairs(i,2) = find((t<=OFF_periods(i,2)),1,'last')+OffSet; % find the end poingt for OFF period i
end
OFF_periods(1,1)=OFF_periods(1,2)-0.5;
OFF_avi_pairs(1,1) = find((t>=OFF_periods(1,1)),1,'first')+OffSet; % find the start point for OFF period 1
OFF_tdms_pairs(1,1) = round(OFF_tdms_pairs(1,2)-(OFF_periods(1,2)-OFF_periods(1,1))/OFF_periods(1,2)*OFF_tdms_pairs(1,2)); % find the start point for OFF period 1

% set rec pairs for avi
REC_avi_pairs=zeros(size(REC_tdms_pairs,1),2); % OFF pairs in the video 
for i=1: size(REC_tdms_pairs,1)
REC_avi_pairs(i,1) = find((t>=REC_periods(i,1)),1,'first')+OffSet; % find the start point for REC period i
REC_avi_pairs(i,2) = find((t<=REC_periods(i,2)),1,'last')+OffSet; % find the end poingt for REC period i
end

%% Combine
len_avi=length(freq_avi);
len_tdms=length(freq_tdms);
if len_avi>len_tdms
    freq_avi(len_tdms+1:len_avi)=[];
    coepow_avi(len_tdms+1:len_avi,:)=[];
else
    freq_tdms(len_avi+1:len_tdms)=[];
    coepow_tdms(len_avi+1:len_tdms,:)=[];
end
%% select dissociation and recovery from WT
%off
coepow_tdms_off=[]; % store the coefficients of WT for dissociation 
coepow_avi_off=[];
realtime_tdms_off=[]; % correspond the real time for each period
realtime_avi_off=[];
inter_tdms_off=[]; % dashed lines inserted to denote each period
inter_avi_off=[];
wt_avi_step=1; % don't change this
wt_tdms_step=tdms_subsamp; %don;t change this
for i=1: size(OFF_periods,1)
coepow_tdms_off=[coepow_tdms_off ...
    coepow_tdms(:,ceil(OFF_tdms_pairs(i,1)/wt_tdms_step):ceil(OFF_tdms_pairs(i,2)/wt_tdms_step))]; 
coepow_avi_off=[coepow_avi_off ...
    coepow_avi(:,ceil(OFF_avi_pairs(i,1)/wt_avi_step):ceil(OFF_avi_pairs(i,2)/wt_avi_step))]; 
% tt=linspace(OFF_periods(i,1),OFF_periods(i,2),OFF_tdms_pairs(i,2)-OFF_tdms_pairs(i,1)+1); 
% realtime_tdms_off=[realtime_tdms_off tt]; 
% inter_tdms_off=[inter_tdms_off length(realtime_tdms_off)]; 
tt=linspace(OFF_periods(i,1),OFF_periods(i,2),OFF_avi_pairs(i,2)-OFF_avi_pairs(i,1)+1); 
realtime_avi_off=[realtime_avi_off tt]; 
inter_avi_off=[inter_avi_off length(realtime_avi_off)]; 
end;

%rec
coepow_tdms_rec=[]; % store the coefficients of WT for dissociation 
coepow_avi_rec=[];
realtime_tdms_rec=[]; % correspond the real time for each period
realtime_avi_rec=[];
inter_tdms_rec=[]; % dashed lines inserted to denote each period
inter_avi_rec=[];
for i=1: size(REC_periods,1)
coepow_tdms_rec=[coepow_tdms_rec ...
    coepow_tdms(:,ceil(REC_tdms_pairs(i,1)/wt_tdms_step):ceil(REC_tdms_pairs(i,2)/wt_tdms_step))]; 
coepow_avi_rec=[coepow_avi_rec ...
    coepow_avi(:,ceil(REC_avi_pairs(i,1)/wt_avi_step):ceil(REC_avi_pairs(i,2)/wt_avi_step))]; 
% tt=linspace(REC_periods(i,1),REC_periods(i,2),REC_tdms_pairs(i,2)-REC_tdms_pairs(i,1)+1); 
% realtime_tdms_rec=[realtime_tdms_rec tt]; 
% inter_tdms_rec=[inter_tdms_rec length(realtime_tdms_rec)]; 
tt=linspace(REC_periods(i,1),REC_periods(i,2),REC_avi_pairs(i,2)-REC_avi_pairs(i,1)+1); 
realtime_avi_rec=[realtime_avi_rec tt]; 
inter_avi_rec=[inter_avi_rec length(realtime_avi_rec)]; 
end;

%% Combine 

co_freq=(freq_tdms+freq_avi)./2;
% off 
co_tdms_off=zeros(size(coepow_avi_off,1),size(coepow_avi_off,2));
div=size(coepow_tdms_off,2)/size(coepow_avi_off,2);
for i =1: size(co_tdms_off,2) 
    co_tdms_off(:,i)=mean(coepow_tdms_off(:,(i-1)*div+1:(i-1)*div+div),2);
end
% normalization
co_tdms_off=co_tdms_off./max(max(co_tdms_off));
coepow_avi_off=coepow_avi_off./max(max(coepow_avi_off));
co_off=coepow_avi_off+co_tdms_off;
co_off=co_off./max(max(co_off));
% rec 
co_tdms_rec=zeros(size(coepow_avi_rec,1),size(coepow_avi_rec,2));
div=size(coepow_tdms_rec,2)/size(coepow_avi_rec,2);
for i =1: size(co_tdms_rec,2) 
    co_tdms_rec(:,i)=mean(coepow_tdms_rec(:,(i-1)*div+1:(i-1)*div+div),2);
end
% normalization
co_tdms_rec=co_tdms_rec./max(max(co_tdms_rec));
coepow_avi_rec=coepow_avi_rec./max(max(coepow_avi_rec));
co_rec=coepow_avi_rec+co_tdms_rec;
co_rec=co_rec./max(max(co_rec));
rt_rec=round(REC_periods-OFF_periods(1,1)+round(OFF_cnt/fs)/2);
rt_off=round(OFF_periods-OFF_periods(1,1)+round(OFF_cnt/fs)/2);
file_out = strcat(fileName(1:end-5),'-COMWAV.mat');
save(file_out,'coepow_avi_rec','co_tdms_rec','co_rec','coepow_avi_off',...
    'co_tdms_off','co_tdms_rec','co_off','inter_avi_off','inter_avi_rec','realtime_avi_off','realtime_avi_rec',...
    'freq_tdms','freq_avi','co_freq','rt_rec','rt_off','wt_tdms_step','REC_periods','REC_array','REC_avi_pairs',...
    'REC_tdms_pairs','OffSet','OFF_tdms_pairs','OFF_periods','OFF_avi_pairs','OFF_cnt','fs');

