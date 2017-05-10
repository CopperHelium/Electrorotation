clc 
close all;
clear all;


% select '*-PeakSteps.mat'
start_path=pwd;
data_chose = uigetfile(start_path,'Select the *-PeakSteps.mat');
data_path=strcat(start_path,'\',data_chose);
load(data_path);

%% Finding peaks for off 
y=real_peak_off;
len=length(real_peak_off);
x=1:1:len;
%x=realtime_avi_off;
%Lower bounds
xl=1;
%Upper bounds
xu=len;
%store the location of step
step=[];
%the size of the step
step_height=[];

% 
threshold=1e3;
stepIteration=5;
Nw=threshold+1;
i=0;
% while max(Nw)>threshold
for ii=1:stepIteration
intnum=length(xl);
step_temp=[];
for i=1:intnum
    step1=stepfit(x(xl(i):xu(i)),y(xl(i):xu(i)));
    step_temp=[step_temp step1'];
end

Nw=xu-xl+1;%the step ‘window’
Delta=abs(step_temp(1,:)); %fitted step size
Q=Delta.*(Nw).^0.5; %parameter to help us to choose best fitted step
[q k]=max(Q);
%get the three parameters for the fitting funtion,
%The first row is step size, the second row is step location.
step=[step step_temp(:,k)]; 
[steppos r]=sort(step(2,:));
step=step(:,r);

%divide into several plateaus
xxl=ceil(step(2,:));
xxu=floor(step(2,:));

xl=sort([1 xxl]);
xu=sort([len xxu]);

if length(step_height)<1
    step_height=[step_temp(3,k)-step_temp(1,k) step_temp(3,k)+step_temp(1,k)];
else
if (k>1)&&(k<length(step_height))
    step_height=[step_height(1:k-1) step_temp(3,k)-step_temp(1,k) step_temp(3,k)+step_temp(1,k) step_height(k+1:end)];
elseif k==1
    step_height=[step_temp(3,k)-step_temp(1,k) step_temp(3,k)+step_temp(1,k) step_height(k+1:end)];
else
    step_height=[step_height(1:k-1) step_temp(3,k)-step_temp(1,k) step_temp(3,k)+step_temp(1,k)];
end
end

end

xstep=sort([xl xu]);
xstep_height=zeros(1, 2*length(step_height));
xstep_height(1:2:2*length(step_height)-1)=step_height;
xstep_height(2:2:2*length(step_height))=step_height;


hold on;
plot(x,y,'o--','color',[0.8 0.8 0.8]);
label_step=ceil(len/6);% the interbal of x axis label for off periods.
set(gca,'xtick',[]);
set(gca,'xtick',x(1:label_step:end));
set(gca,'xticklabel',{round(realtime_avi_off(1:label_step:end))});  

plot(xstep,xstep_height,'color','k','LineWidth',2);
%plot(xstep,xstep_height,'color','k','LineWidth',2);
xlabel('Time (s)');
ylabel('Frequency (Hz)');

xlim([min(x) max(x)]);
box on;


%% Finding peaks for REC
y=real_peak_rec;
len=length(real_peak_rec);
x=1:1:len;
%x=realtime_avi_off;
%Lower bounds
xl=1;
%Upper bounds
xu=len;
%store the location of step
step=[];
%the size of the step
step_height=[];

% 
threshold=1e3;
stepIteration=5;
Nw=threshold+1;
i=0;
% while max(Nw)>threshold
for ii=1:stepIteration
intnum=length(xl);
step_temp=[];
for i=1:intnum
    step1=stepfit(x(xl(i):xu(i)),y(xl(i):xu(i)));
    step_temp=[step_temp step1'];
end

Nw=xu-xl+1;%the step ‘window’
Delta=abs(step_temp(1,:)); %fitted step size
Q=Delta.*(Nw).^0.5; %parameter to help us to choose best fitted step
[q k]=max(Q);
%get the three parameters for the fitting funtion,
%The first row is step size, the second row is step location.
step=[step step_temp(:,k)]; 
[steppos r]=sort(step(2,:));
step=step(:,r);

%divide into several plateaus
xxl=ceil(step(2,:));
xxu=floor(step(2,:));

xl=sort([1 xxl]);
xu=sort([len xxu]);

if length(step_height)<1
    step_height=[step_temp(3,k)-step_temp(1,k) step_temp(3,k)+step_temp(1,k)];
else
if (k>1)&&(k<length(step_height))
    step_height=[step_height(1:k-1) step_temp(3,k)-step_temp(1,k) step_temp(3,k)+step_temp(1,k) step_height(k+1:end)];
elseif k==1
    step_height=[step_temp(3,k)-step_temp(1,k) step_temp(3,k)+step_temp(1,k) step_height(k+1:end)];
else
    step_height=[step_height(1:k-1) step_temp(3,k)-step_temp(1,k) step_temp(3,k)+step_temp(1,k)];
end
end

end

xstep=sort([xl xu]);
xstep_height=zeros(1, 2*length(step_height));
xstep_height(1:2:2*length(step_height)-1)=step_height;
xstep_height(2:2:2*length(step_height))=step_height;


hold on;
plot(x,y,'o--','color',[0.8 0.8 0.8]);
label_step=ceil(len/6);% the interbal of x axis label for off periods.
set(gca,'xtick',[]);
set(gca,'xtick',x(1:label_step:end));
set(gca,'xticklabel',{round(realtime_avi_rec(1:label_step:end))});  

plot(xstep,xstep_height,'color','k','LineWidth',2);
%plot(xstep,xstep_height,'color','k','LineWidth',2);
xlabel('Time (s)');
ylabel('Frequency (Hz)');

xlim([min(x) max(x)]);
box on;

