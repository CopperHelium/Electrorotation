clc 
close all;
clear all;


% select '*-COMWAV.mat'
start_path=pwd;
data_chose = uigetfile(start_path,'Select the *-COMWAV.mat');
data_path=strcat(start_path,'\',data_chose);
load(data_path);
load([data_chose(1:end-10) 'PeakSteps.mat'])

h = figure;
    A1 = axes('position',[0.08 0.2 0.4 0.6]);
    A2 = axes('position',[0.56 0.2 0.4 0.6]);

% Dissociation 
spec_sampling=1; 

set(gcf,'CurrentAxes',A1)
x=(linspace(1,size(co_off(:,1:spec_sampling:end),2),size(co_off(:,1:spec_sampling:end),2)));
y=co_freq(1:end);
[X Y]=meshgrid(x,y);

   % plot wavelet spectrum
   surf(X, Y,co_off(:,1:spec_sampling:end),'LineStyle','none');
   c=colorbar;
   c.Label.String = 'Intensity';
   clim_min=0.1;
   clim_max=1;
   set(gca,'Clim',[clim_min clim_max]);
   %plot period lines
  xlabel('Time (s)');
  ylabel('Frequency (Hz)');
  plot_name= data_chose(1:end-11);
  plot_name(find(plot_name==char('_')))=[];
  name=['Combined Dissociation' ',' plot_name ',' 'Iteration:' '1s'];
  title(name);
  axis tight;
 % set x label for each pr
  label_step=ceil(size(OFF_avi_pairs,1)/6); % the interbal of x axis label for off periods.
   % plot peak  
  inter_avi_off2=ceil(inter_avi_off./spec_sampling);
  z=40*ones(1,length(real_peak_off));
  hold on;
  peaknum = linspace(1,length(real_peak_off),length(real_peak_off));
  p1=plot3(peaknum(1:end),real_peak_off(1:end),z(1:end),'.','color','k','MarkerSize',8);
  %p2=plot3(peaknum(1:end),pottsL2_off,z,'g--.', 'MarkerSize', 10);
  p3=plot3(peaknum(1:end),step_modified_off,z,'r--.', 'MarkerSize', 10);
  legend([p1 p3],'Peaks','Steps','Location','best');
  set(gca,'xtick',[]);
  set(gca,'xtick',inter_avi_off2(2:label_step:end)-floor((inter_avi_off2(2)-inter_avi_off2(1))/2)+1);
  shift_time=OFF_periods(1,1);
  set(gca,'xticklabel',{round(OFF_periods(2:label_step:end,1)-shift_time+round(OFF_cnt/fs)/2)});  
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
spec_sampling=1; % sampling frequency for the wavelet spectrum
x=(linspace(1,size(co_rec(:,1:spec_sampling:end),2),size(co_rec(:,1:spec_sampling:end),2)));
y=co_freq(1:end);
[X Y]=meshgrid(x,y);


   % plot wavelet spectrum
   surf(X, Y,co_rec(:,1:spec_sampling:end),'LineStyle','none');
   c=colorbar;
   c.Label.String = 'Intensity';
   clim_min=0.1;
   clim_max=1;
   set(gca,'Clim',[clim_min clim_max]);
   
  xlabel('Time (s)');
  ylabel('Frequency (Hz)');
  name=['Combined Recovery' ',' plot_name];
  title(name);
  axis tight;
 % set x label for each pr

 % plot peak  
  z=40*ones(1,length(real_peak_rec));
  inter_avi_rec2=floor(inter_avi_rec./spec_sampling);
  hold on;
  peaknum = linspace(1,length(real_peak_rec),length(real_peak_rec));
  p1=plot3(peaknum(1:end),real_peak_rec(1:end),z(1:end),'.','color','k','MarkerSize',8);
  %p2=plot3(peaknum(1:end),pottsL2_rec,z,'g--.', 'MarkerSize', 10);
  p3=plot3(peaknum(1:end),step_modified_rec,z,'r--.', 'MarkerSize', 10);
  view([0 90]);
  legend([p1 p3],'Peaks','Steps','Location','best');
  label_step=ceil(size(REC_periods,1)/6);% the interbal of x axis label for off periods.
  set(gca,'xtick',[]);
  set(gca,'xtick',inter_avi_rec2(1:label_step:end)-floor((inter_avi_rec2(2)-inter_avi_rec2(1))/2)+1);
  set(gca,'xticklabel',{round(REC_periods(1:label_step:end,1)-shift_time+round(OFF_cnt/fs)/2)});  

% save png file  
set(gcf,'PaperType','A4')
    set(gcf,'PaperUnits','centimeters')
    xSize = 28; ySize = 15;
    xLeft = (30-xSize)/2; yTop = (21-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[100 100 xSize*50 ySize*50])
    set(gcf,'PaperOrientation','landscape');
    plot_name = strcat(tdms_name(1:end-4),'combine.png');
    saveas(gcf,plot_name)
