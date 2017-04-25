% This code is for combination of tdms and avi files 
clc 
close all;
clear all;
%% load tdms and avi files 
tdms_name= '16h_30m_14s.mat';
video_path='6.avi';
avi_name=[tdms_name(1:end-4) '_video.mat'];
load(tdms_name);
load(avi_name);
obj = VideoReader(video_path);
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
    else
        OFF_tdms_pairs = vertcat([1 OnOff(1)],OFF_tdms_pairs);
    end
     REC_array = OFF_tdms_pairs(end,2):OFF_cnt:size(X_signal,1);
     REC_tdms_pairs = [REC_array(1:end-1)' REC_array(2:end)'];
     OFF_periods= (OFF_tdms_pairs-1)./fs;
     REC_periods=(REC_tdms_pairs-1)./fs;

%% WT for avi file, only preserve interested area
vidHeight = obj.Height;
vidWidth = obj.Width;
vidFrames = obj.Duration*obj.FrameRate;
T = obj.Duration/vidFrames;           
Fs = obj.FrameRate;   
t = (0:vidFrames-1)*T;        
wlet = 'morse';              
Scale_avi_min=0;  
Scale_avi_max=1e7; 
wt_avi_step=1;
[coefsx, freqx] = cwt(pos(1:wt_avi_step:end,1),wlet,Fs,'TimeBandwidth',60,'VoicesPerOctave',48); 
[coefsy, freqy] = cwt(pos(1:wt_avi_step:end,2),wlet,Fs,'TimeBandwidth',60,'VoicesPerOctave',48);
Scale_avi_min=max(Scale_avi_min,find(freqx>20, 1, 'last' ));
Scale_avi_max=min(Scale_avi_max,length(freqx));
coepow_avi=abs(coefsy)+abs(coefsx);
coepow_avi(Scale_avi_max+1:end,:)=[];
coepow_avi(1:Scale_avi_min-1,:)=[];
freq_avi =freqx(Scale_avi_min:1:Scale_avi_max);
clear coefsx;
clear coefsy;
% coepow_avi = coepow;
% freq_avi=freq;
% clear coepow;
% clear freq;
%% Find OFF and REC pairs by related tdms files
% determine offset first
value_min=0.04; 
pre_peak=wltpeakfinder(coepow_avi(:,1:20*Fs),3,value_min);
pre_peak(isnan(pre_peak))=0;
pre_switch=diff(pre_peak);
off_start=find(abs(pre_switch)>30,1,'first');
OffSet=off_start-floor(OFF_periods(1,2)*Fs)-1;
% set off pairs for avi
OFF_avi_pairs=zeros(size(OFF_tdms_pairs,1),2); % OFF pairs in the video 
for i=1: size(OFF_tdms_pairs,1)
OFF_avi_pairs(i,1) = find((t>=OFF_periods(i,1)),1,'first')+OffSet; % find the start point for OFF period i
OFF_avi_pairs(i,2) = find((t<=OFF_periods(i,2)),1,'last')+OffSet; % find the end poingt for OFF period i
end
% set rec pairs for avi
REC_avi_pairs=zeros(size(REC_tdms_pairs,1),2); % OFF pairs in the video 
for i=1: size(REC_tdms_pairs,1)
REC_avi_pairs(i,1) = find((t>=REC_periods(i,1)),1,'first')+OffSet; % find the start point for REC period i
REC_avi_pairs(i,2) = find((t<=REC_periods(i,2)),1,'last')+OffSet; % find the end poingt for REC period i
end

%% WT for tdms file, only preserve interested area
Scale_tdms_min=0;
Scale_tdms_max=1e7;
wt_tdms_step=20;
wlet='morse';
[coefsx, freqx] = cwt(X_signal(1:wt_tdms_step:end),wlet,fs/wt_tdms_step,'TimeBandwidth',60,'VoicesPerOctave',48); 
[coefsy, freqy] = cwt(Y_signal(1:wt_tdms_step:end),wlet,fs/wt_tdms_step,'TimeBandwidth',60,'VoicesPerOctave',48);
Scale_tdms_min=max(Scale_tdms_min,find(freqx>20, 1, 'last' ));
Scale_tdms_max=min(Scale_tdms_max,length(freqx));
coepow_tdms=abs(coefsy)+abs(coefsx);
coepow_tdms(Scale_tdms_max+1:end,:)=[];
coepow_tdms(1:Scale_tdms_min-1,:)=[];
freq_tmds =freqx(Scale_tdms_min:1:Scale_tdms_max);
clear coefsx;
clear coefsy;

len_avi=length(freq_avi);
len_tdms=length(freq_tmds);
if len_avi>len_tdms
    freq_avi(len_tdms+1:len_avi)=[];
    coepow_avi(len_tdms+1:len_avi,:)=[];
else
    freq_tmds(len_avi+1:len_tmds)=[];
    coepow_tmds(len_avi+1:len_tmds,:)=[];
end
%% select dissociation and recovery from WT
%off
coepow_tdms_off=[]; % store the coefficients of WT for dissociation 
coepow_avi_off=[];
realtime_tdms_off=[]; % correspond the real time for each period
realtime_avi_off=[];
inter_tdms_off=[]; % dashed lines inserted to denote each period
inter_avi_off=[];
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
inter_avi_off=inter_avi_off-1;
inter_avi_off=[0 inter_avi_off];
inter_avi_off(end)=[];
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

co_freq=(freq_tmds+freq_avi)./2;
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

%% Peak finder and step finder 
% OFF
        value_min=0.5; % threshold of peak
        peak_sampling_off= 1; % sampling frequency for peak finder
        real_peak_off=[];
        % a customed function for peak finder
        peak_off=wltpeakfinder(co_off(:,1:peak_sampling_off:end),3,value_min);
       
        % correspond the real value of peaks for the WT scale
        tt = isnan(peak_off);
        for i=1: length(peak_off)
        if tt(i) == 0
            scale= peak_off(i);     
            real_peak_off(i)= co_freq(scale);     
        else  real_peak_off(i)= 0;
        end
        end     
        % minL2Potts: A function uses potts model.
        % Parameter 1: data
        % Parameter 2: level of peak. The bigger number will be more coarse.
        pottsL2_off = minL2Potts(real_peak_off, 10);
        zeropos=find(pottsL2_off==0);
        pottsL2_off(zeropos)=nan;
        % modify the step by "stepmodify" function. The second parameter is the thresh which considered as
        %the same steps;
        step_modified_off = stepmodify(pottsL2_off,0.3,4);
        step_modified_off(find(step_modified_off==0))=nan;
% REC
        value_min=0; % threshold of peak
        peak_sampling_rec= 50; % sampling frequency for peak finder
        real_peak_rec=[];
        % a customed function for peak finder
        peak_rec=wltpeakfinder(co_rec(:,1:peak_sampling_rec:end),3,value_min);
       
        % correspond the real value of peaks for the WT scale
        tt = isnan(peak_rec);
        for i=1: length(peak_rec)
        if tt(i) == 0
            scale= peak_rec(i);     
            real_peak_rec(i)= co_freq(scale);     
        else  real_peak_rec(i)= 0;
        end
        end   
        % minL2Potts: A function uses potts model.
        % Parameter 1: data
        % Parameter 2: level of peak. The bigger number will be more coarse.
        pottsL2_rec = minL2Potts(real_peak_rec, 5);
        zeropos=find(pottsL2_rec==0);
        pottsL2_rec(zeropos)=nan;
        % modify the step by "stepmodify" function. The second parameter is the thresh which considered as
        %the same steps;
        step_modified_rec = stepmodify(pottsL2_rec,0.3,4);
        step_modified_rec(find(step_modified_rec==0))=nan;
        
%% Plot
 h = figure;
    A1 = axes('position',[0.08 0.2 0.4 0.6]);
    A2 = axes('position',[0.56 0.2 0.4 0.6]);

% Dissociation 
% set the grid for the wavelet spectrum
spec_sampling=1; % sampling frequency for the wavelet spectrum

set(gcf,'CurrentAxes',A1)
x=(linspace(1,size(co_off(:,1:spec_sampling:end),2),size(co_off(:,1:spec_sampling:end),2)));
y=co_freq(1:end);
[X Y]=meshgrid(x,y);

   % plot wavelet spectrum
   surf(X, Y,co_off(:,1:spec_sampling:end),'LineStyle','none');
   c=colorbar;
   c.Label.String = 'Intensity';
   clim_min=0.3;
   clim_max=2;
   set(gca,'Clim',[clim_min clim_max]);
   %plot period lines
  xlabel('Time (s)');
  ylabel('Frequency (Hz)');
  tdms_name(find(tdms_name==char('_')))=[];
  name=['Combined Dissociation' ',' tdms_name];
  title(name);
  axis tight;
 % set x label for each pr
  label_step=size(OFF_avi_pairs,1)/5; % the interbal of x axis label for off periods.
   % plot peak  
  inter_avi_off=ceil(inter_avi_off./spec_sampling);
  z=40*ones(1,length(peak_off));
  hold on;
  peaknum = linspace(1,length(peak_off),length(peak_off));
  plot3(peaknum(1:end),real_peak_off(1:end),z(1:end),'.','color','k','MarkerSize',8);
  plot3(peaknum(1:end),pottsL2_off,z,'g--.', 'MarkerSize', 10);
  plot3(peaknum(1:end),step_modified_off,z,'r--.', 'MarkerSize', 10);
  legend('Intensity of WT','Peak','Step','Step modified','Location','best');
  set(gca,'xtick',[]);
  set(gca,'xtick',inter_avi_off(1:label_step:end));
  set(gca,'xticklabel',{floor(OFF_periods(1:label_step:end,1))});  
   yline=[min(co_freq) max(co_freq)];
   zline=[20 20];
   xline=[0 0];
   hold on; 
for i=1:length(inter_avi_off)
    plot3(xline+inter_avi_off(i),yline,zline,'--w');
end
  view([0 90]);
  
  
% Recovery
% set the grid for the wavelet spectrum
set(gcf,'CurrentAxes',A2)
spec_sampling=peak_sampling_rec; % sampling frequency for the wavelet spectrum
x=(linspace(1,size(co_rec(:,1:spec_sampling:end),2),size(co_rec(:,1:spec_sampling:end),2)));
y=co_freq(1:end);
[X Y]=meshgrid(x,y);


   % plot wavelet spectrum
   surf(X, Y,co_rec(:,1:spec_sampling:end),'LineStyle','none');
   c=colorbar;
   c.Label.String = 'Intensity';
   clim_min=0.3;
   clim_max=2;
   set(gca,'Clim',[clim_min clim_max]);
   
  xlabel('Time (s)');
  ylabel('Frequency (Hz)');
  name=['Combined Recovery' ',' tdms_name];
  title(name);
  axis tight;
 % set x label for each pr

 % plot peak  
  z=40*ones(1,length(peak_rec));
  inter_avi_rec=floor(inter_avi_rec./spec_sampling);
  hold on;
  peaknum = linspace(1,length(peak_rec),length(peak_rec));
  plot3(peaknum(1:end),real_peak_rec(1:end),z(1:end),'.','color','k','MarkerSize',8);
  plot3(peaknum(1:end),pottsL2_rec,z,'g--.', 'MarkerSize', 10);
  plot3(peaknum(1:end),step_modified_rec,z,'r--.', 'MarkerSize', 10);
  view([0 90]);
  legend('Intensity of WT','Peak','Step','Step modified','Location','best');
  label_step=size(REC_periods,1)/6;% the interbal of x axis label for off periods.
  set(gca,'xtick',[]);
  set(gca,'xtick',inter_avi_rec(1:label_step:end));
  set(gca,'xticklabel',{floor(REC_periods(1:label_step:end,1))});  

% save png file  
set(gcf,'PaperType','A4')
    set(gcf,'PaperUnits','centimeters')
    xSize = 28; ySize = 15;
    xLeft = (30-xSize)/2; yTop = (21-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[100 100 xSize*50 ySize*50])
    set(gcf,'PaperOrientation','landscape');
    plot_name = strcat(tdms_name(1:end-4),'combine25.png');
    saveas(gcf,plot_name)
