% Run wavelet for the tdms files('*s.mat')

clc 
close all;
clear all;

%start from the current path
start_path=pwd;
%select the '*s.mat' wanted
TDMS_name = uigetfile(start_path,'Select the *s.mat');
TDMS_path=strcat(start_path,'\',TDMS_name);
load(TDMS_path);

 %% Load TDMS and rearrange each REC and OFF pairs.
    X_signal = UntitledVoltage_0.Data; % Contains X-Photomultiplier signal
    Y_signal = UntitledVoltage_1.Data; % Contains Y-Photomultiplier signal
    fs = UntitledVoltage_2.Property.wf_samples; % number of samples per sec
    F_state = UntitledVoltage_2.Data; % Shows whether the field was ON or OFF
    F_state = (F_state > 2); % convert to binary (0/1)
    
    Switch = diff(F_state); % find if the field was switched
    OnOff = find(Switch); % Indices of the switching data points
    
    ON_pairs = [OnOff(1:2:end) OnOff(2:2:end)];
    OFF_pairs = [OnOff(2:2:end-2) OnOff(3:2:end)];
    
    % unify the selection of ON_cnt or OFF_cnt in case of not enough pairs
    numoff=size(OFF_pairs,1);
    if numoff>=3
        numoff=3;
    end
    numon=size(ON_pairs,1);
    if numon>=3
        numon=3;
    end
    ON_cnt = (ON_pairs(numon,2)-ON_pairs(numon,1)); % # samples in each ON event
    OFF_cnt = (OFF_pairs(numoff,2)-OFF_pairs(numoff,1)); % # samples in each OFF event 

    rec_flag = 0;
    if size(F_state,1)-ON_pairs(end,2)>OFF_cnt && F_state(end)==0
        %some data need (OFF_pairs(1,2)-OFF_pairs(1,1))
        OFF_pairs = vertcat([1 OnOff(1)],OFF_pairs,[OnOff(end) OnOff(end)+(OFF_pairs(numoff,2)-OFF_pairs(numoff,1))]);
        rec_flag = 1;
    elseif size(F_state,1)-ON_pairs(end,2)<OFF_cnt && F_state(end)==0
        OFF_pairs = vertcat([1 OnOff(1)],OFF_pairs);
    else
        OFF_pairs = vertcat([1 OnOff(1)],OFF_pairs);
    end
    
     REC_array = OFF_pairs(end,2):OFF_cnt:size(X_signal,1);
     REC_pairs = [REC_array(1:end-1)' REC_array(2:end)'];
     OFF_periods= (OFF_pairs-1)./fs;
     REC_periods=(REC_pairs-1)./fs;
     %% Wavelet Transform part     
% Notice: here for simplification, the sampling have been decreased to
% 1/20;
Scale_min=0;
Scale_max=1e7;
sub_sampling=20;
wlet='morse';
[coefsx, freqx] = cwt(X_signal(1:sub_sampling:end),wlet,fs/sub_sampling,'TimeBandwidth',60,'VoicesPerOctave',48); 
[coefsy, freqy] = cwt(Y_signal(1:sub_sampling:end),wlet,fs/sub_sampling,'TimeBandwidth',60,'VoicesPerOctave',48);
Scale_min=max(Scale_min,find(freqx>20, 1, 'last' ));
Scale_max=min(Scale_max,length(freqx));
coepow=abs(coefsy)+abs(coefsx);
coepow(Scale_max+1:end,:)=[];
coepow(1:Scale_min-1,:)=[];
freq =freqx(Scale_min:1:Scale_max);

%get the time vector
tdmsFrames = length(X_signal);
Time = tdmsFrames./fs;   % the whole time duration
t = (0:sub_sampling:tdmsFrames-1)*(Time./tdmsFrames); % Time vector
  
file_out = strcat(fileName(1:end-5),'-TDMSWAV.mat');
save(file_out,'coepow','freq','t','sub_sampling');

