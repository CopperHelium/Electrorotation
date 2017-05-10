clc 
close all;
clear all;

% select '*-COMWAV.mat'
start_path=pwd;
data_chose = uigetfile(start_path,'Select the *-COMWAV.mat');
data_path=strcat(start_path,'\',data_chose);
load(data_path);



%% Peak finder and step finder 
% OFF
        value_min=0.04; % threshold of peak
        peak_sampling_off= 1; % sampling frequency for peak finder
        real_peak_off=[];
        peak_off=wltpeakfinder(co_off(:,1:peak_sampling_off:end),3,value_min);
       
        tt = isnan(peak_off);
        for i=1: length(peak_off)
        if tt(i) == 0
            scale= peak_off(i);     
            real_peak_off(i)= co_freq(scale);     
        else  real_peak_off(i)= 0;
        end
        end     
      
        pottsL2_off = minL2Potts(real_peak_off, 10);
        zeropos=find(pottsL2_off==0);
        pottsL2_off(zeropos)=nan;
        step_modified_off = stepmodify(pottsL2_off,0.3,8);
        step_modified_off(find(step_modified_off==0))=nan;
% REC
        value_min=0; % threshold of peak
        peak_sampling_rec= 1; % sampling frequency for peak finder
        real_peak_rec=[];
        peak_rec=wltpeakfinder(co_rec(:,1:peak_sampling_rec:end),3,value_min);
        tt = isnan(peak_rec);
        for i=1: length(peak_rec)
        if tt(i) == 0
            scale= peak_rec(i);     
            real_peak_rec(i)= co_freq(scale);     
        else  real_peak_rec(i)= 0;
        end
        end   
        pottsL2_rec = minL2Potts(real_peak_rec, 5);
        zeropos=find(pottsL2_rec==0);
        pottsL2_rec(zeropos)=nan;
        step_modified_rec = stepmodify(pottsL2_rec,0.3,4);
        step_modified_rec(find(step_modified_rec==0))=nan;
        file_out=[data_chose(1:end-10) 'PeakSteps.mat'];
        save(file_out,'inter_avi_off','real_peak_off','inter_avi_rec',...
            'real_peak_rec','realtime_avi_off','realtime_avi_rec',...
            'step_modified_rec','step_modified_off');
