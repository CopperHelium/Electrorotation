% This is a ".m" file to determine the frequency among time domain for tdms
% files. Use the "simpleConvertTDMS.m" or "ConvertTDMS.m" first to convert
% the TDMS files to mat format. To determine the steps, pottslab package is
% used here so JavaPath should be set first by run the "installPottslab.m"
% in the "Pottslab0.5" folder.This version include analysis of dissociation
% and recovery periods.

clc;
close all;
clear all;
%% Load Data by setting path and selecting data folder
warning('off', 'MATLAB:colon:nonIntegerIndex')

% Ask the user to select folder containing the .mat data files to be
% analyzed for freq. Use current working directory as a starting point.

start_path=pwd;
folder_name = uigetdir(start_path,'Select folder containing .mat files to be analyzed:');
cd(folder_name)
%----------------------------------------------------------------

%Collect a list of all the files in the current directory
Files=dir('*s.mat');

%Measure how long the list of files is
lengthFiles=length(Files);


%% Run for each data file separately

   for kk = 1:lengthFiles
    %clearvars -except kk Files lengthFiles folder_name
    clc;
 
    
    file_name = Files(kk).name;
    path_name = strcat(folder_name,'\',file_name);
    
    %delete "_" in the file name.
    
    
    %Load the mat file containing the data. Eventually this should be in a loop
    %going through all files in a cerain folder(s).
    load(path_name)
    fileName = file_name;
    fileName(find(fileName==char('_')))=[];
    
   %% --------------%% Start the main codes %%---------------------------------------
    %% Load measured parameters and arrange them to pairs.
    
    X_signal = UntitledVoltage_0.Data; % Contains X-Photomultiplier signal
    Y_signal = UntitledVoltage_1.Data; % Contains Y-Photomultiplier signal
    fs = UntitledVoltage_2.Property.wf_samples; % number of samples per sec
    F_state = UntitledVoltage_2.Data; % Shows whether the field was ON or OFF
    F_state = (F_state > 2); % convert to binary (0/1)
    
    Switch = diff(F_state); % find if the field was switched
    OnOff = find(Switch); % Indices of the switching data points
    
    %Define pairs in sample numbers describing the times when the field was on
    %and off. The first column is when the event startes and the second when it
    %ended. Frequency calculated between these two.
    
    ON_pairs = [OnOff(1:2:end) OnOff(2:2:end)];
    OFF_pairs = [OnOff(2:2:end-2) OnOff(3:2:end)];
    
    
    numoff=size(OFF_pairs,1);
    if numoff>=3
        numoff=3;
    end
    numon=size(ON_pairs,1);
    if numon>=3
        numon=3;
    end
    ON_cnt = (ON_pairs(numon,2)-ON_pairs(numon,1)); % # samples in each ON event
    
    %some data need (OFF_pairs(1,2)-OFF_pairs(1,1))
    OFF_cnt = (OFF_pairs(numoff,2)-OFF_pairs(numoff,1)); % # samples in each OFF event 
    

    % The matrix above missed the last OFF event immediately following the last
    % ON event and before the recovery starts (two seconds after the last ON).
    % Adding that bit to the Off_pairs here.
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
     %% WT part
Scale_min=0;
Scale_max=1e7;
wave_step=10;
wlet='morse';
[coefsx, freqx] = cwt(X_signal(1:wave_step:end),wlet,fs/wave_step,'TimeBandwidth',60,'VoicesPerOctave',48); 
[coefsy, freqy] = cwt(Y_signal(1:wave_step:end),wlet,fs/wave_step,'TimeBandwidth',60,'VoicesPerOctave',48);
Scale_min=max(Scale_min,find(freqx>20, 1, 'last' ));
Scale_max=min(Scale_max,length(freqx));
coepow=abs(coefsy)+abs(coefsx);
coepow(Scale_max+1:end,:)=[];
coepow(1:Scale_min-1,:)=[];
freq =freqx(Scale_min:1:Scale_max);
   
    %% ------------- %% Dissociation Periods %% --------------------------
%% Dealing with dissociation period

coepow_off=[]; % store the coefficients of WT for dissociation 
real_time_off=[]; % correspond the real time for each period
inter_off=[]; % dashed lines inserted to denote each period

for i=1: size(OFF_pairs,1)
coepow_off=[coepow_off coepow(:,ceil(OFF_pairs(i,1)/wave_step):ceil(OFF_pairs(i,2)/wave_step))]; % get the coefficients for the OFF period i
tt=linspace(OFF_periods(i,1),OFF_periods(i,2),OFF_pairs(i,2)-OFF_pairs(i,1)+1); % get the real time
real_time_off=[real_time_off tt]; % collect real time for fields off
inter_off=[inter_off length(real_time_off)]; % collect lines position 
end;
 inter_off=inter_off-1;
 inter_off=[0 inter_off]; % insert the start time.
 inter_off(end)=[];     
 %% peak finder for OFF
        value_min=1; % threshold of peak
        peak_sampling_off= 1; % sampling frequency for peak finder
        real_peak_off=[];
        % a customed function for peak finder
        peak_off=wltpeakfinder(coepow_off(:,1:peak_sampling_off:end),3,value_min);
       
        % correspond the real value of peaks for the WT scale
        tt = isnan(peak_off);
        for i=1: length(peak_off)
        if tt(i) == 0
            scale= peak_off(i);     
            real_peak_off(i)= freq(scale);     
        else  real_peak_off(i)= 0;
        end
        end  
 %% Step finder for off      
        % minL2Potts: A function uses potts model.
        % Parameter 1: data
        % Parameter 2: level of peak. The bigger number will be more coarse.
        pottsL2_off = minL2Potts(real_peak_off, 20);
        zeropos=find(pottsL2_off==0);
        pottsL2_off(zeropos)=nan;
        % modify the step by "stepmodify" function. The second parameter is the thresh which considered as
        %the same steps;
        step_modified_off = stepmodify(pottsL2_off,0.3,4);
        step_modified_off(find(step_modified_off==0))=nan;

    %% ------------- %% Recovery Periods  %% ----------------------------     
    coepow_rec=[]; % store the coefficients of WT for recovery 
    real_time_rec=[]; % correspond the real time for each period
    inter_rec=[];
    for i=1: size(REC_pairs,1)
    coepow_rec=[coepow_rec coepow(:,ceil(REC_pairs(i,1)./wave_step):ceil(REC_pairs(i,2)./wave_step))]; % get the coefficients for the REC period i
    tt=linspace(REC_periods(i,1),REC_periods(i,2),REC_pairs(i,2)-REC_pairs(i,1)+1); % get the real time
    real_time_rec=[real_time_rec tt]; % collect real time for fields recovery
    inter_rec=[inter_rec length(real_time_rec)]; % collect lines position 
    end;
%% peak finder for REC
        value_min=2; % threshold of peak
        peak_sampling_rec= 50; % sampling frequency for peak finder
        real_peak_rec=[];
        % a customed function for peak finder
        peak_rec=wltpeakfinder(coepow_rec(:,1:peak_sampling_rec:end),3,value_min);
       
        % correspond the real value of peaks for the WT scale
        tt = isnan(peak_rec);
        for i=1: length(peak_rec)
        if tt(i) == 0
            scale= peak_rec(i);     
            real_peak_rec(i)= freq(scale);     
        else  real_peak_rec(i)= 0;
        end
        end
  %% Step finder for off      
        % minL2Potts: A function uses potts model.
        % Parameter 1: data
        % Parameter 2: level of peak. The bigger number will be more coarse.
        pottsL2_rec = minL2Potts(real_peak_rec, 20);
        zeropos=find(pottsL2_rec==0);
        pottsL2_rec(zeropos)=nan;
        % modify the step by "stepmodify" function. The second parameter is the thresh which considered as
        %the same steps;
        step_modified_rec = stepmodify(pottsL2_rec,0.3,4);
        step_modified_rec(find(step_modified_rec==0))=nan;
   %% Plot the dissociation and recovery parts together
    h = figure;
    A1 = axes('position',[0.08 0.2 0.4 0.6]);
    A2 = axes('position',[0.56 0.2 0.4 0.6]);
    set(gcf,'CurrentAxes',A1)
    spec_sampling=1;
    wavelet_spectum(coepow_off(:,1:spec_sampling:end),freq,'Dissociation',fileName);
    ps_plot(wave_step,spec_sampling,OFF_periods,peak_off,real_peak_off,pottsL2_off,step_modified_off,freq,inter_off,real_time_off,size(OFF_pairs,1)/2,0,14,fileName,0);
    view([0 90]);
    
    spec_sampling=50;
    set(gcf,'CurrentAxes',A2)
    wavelet_spectum(coepow_rec(:,1:spec_sampling:end),freq,'Recovery',fileName);
    ps_plot(wave_step,spec_sampling,REC_periods,peak_rec,real_peak_rec,pottsL2_rec,step_modified_rec,freq,inter_rec,real_time_rec,size(REC_pairs,1)/4,0,14,fileName,1);
    view([0 90]);
    
    set(gcf,'PaperType','A4')
    set(gcf,'PaperUnits','centimeters')
    xSize = 28; ySize = 15;
    xLeft = (30-xSize)/2; yTop = (21-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[100 100 xSize*50 ySize*50])
    set(gcf,'PaperOrientation','landscape');
    plot_name = strcat(file_name(1:end-5),'date21.png');
    saveas(gcf,plot_name)
  
   end;


