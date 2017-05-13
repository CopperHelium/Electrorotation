% This is the matlab code to use chi-squared step finding method. Here I
% use non-linear least-sqaured curve fitting to achieve the Chi-squared
% minimization. The fitting type function is actan instead of steo which is
% robust with the initial value. The iteration times are depended on the
% threshold set by hand. 

clc 
close all;
clear all;


% select folder and '*-PeakSteps.mat' file, set the path of data and load
start_path=pwd;
folder_chose= uigetdir(start_path,'Select the folder wanted');
data_chose = uigetfile(folder_chose,'Select the *-PeakSteps.mat');
data_path=strcat(folder_chose,'\',data_chose);
load(data_path);

%% Finding peaks for off 
Measured_Y=real_peak_off;
len=length(real_peak_off);
NumVector=1:1:len;          %just use the simple numbers of the frame.

Step_Left=1;   % The left position of the plateau of step pairs.
Step_Right=len; % The right position of the plateau of step pairs.

Step_Var=[];  % store the parameters of step fitting
Step_Height=[]; % store the height of each plateau 


threshold=3; % Set the threshold to 3. 
Q=1e3; % Set the initial value for quality of the fitting

while max(Q)>threshold  % While the maximum of quality value larger than the threshold, do the loop.
iternum=length(Step_Left); % The num of iteration
step_temp=[]; % store the temporary parameter for step fitting.
for i=1:iternum  % Finding step for each plateau 
    
    % a customed function for least-square fitting 
    step1=stepfit_test(NumVector(Step_Left(i):Step_Right(i)),Measured_Y(Step_Left(i):Step_Right(i)),len);
    % store the parameter of the all parameters of fitting
    step_temp=[step_temp step1'];
end

Nw=Step_Right-Step_Left+1; % the window of each step as explained in the literature.
Delta=abs(step_temp(1,:)); % fitted step size
Q=Delta.*(Nw).^(1/3); % quality value to help us to choose best fitted step
[q k]=max(Q); % get the maximum of Q and the number of fitted step

Step_Var=[Step_Var step_temp(:,k)]; % store the selected fitting parameters

% rearrange the parameter by the the sequence of location of step. 
[steppos r]=sort(Step_Var(2,:));
Step_Var=Step_Var(:,r);

%Divide into separated plateaus
xxl=ceil(Step_Var(2,:));
xxu=floor(Step_Var(2,:));

% rearrange all the plateaus
Step_Left=sort([1 xxl]);
Step_Right=sort([len xxu]);

% Store the step height
if length(Step_Height)<1
    Step_Height=[step_temp(3,k)+(-pi/2)*step_temp(1,k) step_temp(3,k)+pi/2*step_temp(1,k)];
else
if (k>1)&&(k<length(Step_Height))
    Step_Height=[Step_Height(1:k-1) step_temp(3,k)+(-pi/2)*step_temp(1,k) step_temp(3,k)+pi/2*step_temp(1,k) Step_Height(k+1:end)];
elseif k==1
    Step_Height=[step_temp(3,k)+(-pi/2)*step_temp(1,k) step_temp(3,k)+pi/2*step_temp(1,k) Step_Height(k+1:end)];
else
    Step_Height=[Step_Height(1:k-1) step_temp(3,k)+(-pi/2)*step_temp(1,k) step_temp(3,k)+pi/2*step_temp(1,k)];
end
end

end

% sort all the step and set teh height for the steps
xstep=sort([Step_Left Step_Right]);
xstep_height=zeros(1, 2*length(Step_Height));
xstep_height(1:2:2*length(Step_Height)-1)=Step_Height;
xstep_height(2:2:2*length(Step_Height))=Step_Height;


% plot figure.
hold on;
plot(NumVector,Measured_Y,'o--','color',[0.8 0.8 0.8]);
label_step=ceil(len/6);% the interbal of x axis label for off periods.
set(gca,'xtick',[]);
set(gca,'xtick',NumVector(1:label_step:end));
set(gca,'xticklabel',{round(realtime_avi_off(1:label_step:end))});  

plot(xstep,xstep_height,'color','k','LineWidth',2);
%plot(xstep,xstep_height,'color','k','LineWidth',2);
xlabel('Time (s)');
ylabel('Frequency (Hz)');

xlim([min(NumVector) max(NumVector)]);
box on;


%% Finding peaks for rec

Measured_Y=real_peak_rec;
len=length(real_peak_rec);
NumVector=1:1:len;    %just use the simple numbers of the frame.

Step_Left=1; % The left position of the plateau of step pairs.
Step_Right=len;  % The right position of the plateau of step pairs.

Step_Var=[]; % store the parameters of step fitting

Step_Height=[];% store the height of each plateau 



threshold=5; % Set the threshold to 5. Smaller number leads to overfitting
Q=1e3; % Set the initial value for quality of the fitting


while max(Q)>threshold
iternum=length(Step_Left);
step_temp=[];
for i=1:iternum
    
    % a customed function for least-square fitting 
    step1=stepfit_test(NumVector(Step_Left(i):Step_Right(i)),Measured_Y(Step_Left(i):Step_Right(i)),len);
    % store the parameter of the all parameters of fitting
    step_temp=[step_temp step1'];
end

Nw=Step_Right-Step_Left+1;% the window of each step as explained in the literature.
Delta=abs(step_temp(1,:)); %fitted step size
Q=Delta.*(Nw).^(1/2); % quality value to help us to choose best fitted step
[q k]=max(Q);  % get the maximum of Q and the number of fitted step

Step_Var=[Step_Var step_temp(:,k)]; 

% rearrange the parameter by the the sequence of location of step. 
[steppos r]=sort(Step_Var(2,:));
Step_Var=Step_Var(:,r);

%divide into several plateaus
xxl=ceil(Step_Var(2,:));
xxu=floor(Step_Var(2,:));

% rearrange all the plateaus
Step_Left=sort([1 xxl]);
Step_Right=sort([len xxu]);

%store the step height
if length(Step_Height)<1
    Step_Height=[step_temp(3,k)+(-pi/2)*step_temp(1,k) step_temp(3,k)+pi/2*step_temp(1,k)];
else
if (k>1)&&(k<length(Step_Height))
    Step_Height=[Step_Height(1:k-1) step_temp(3,k)+(-pi/2)*step_temp(1,k) step_temp(3,k)+pi/2*step_temp(1,k) Step_Height(k+1:end)];
elseif k==1
    Step_Height=[step_temp(3,k)+(-pi/2)*step_temp(1,k) step_temp(3,k)+pi/2*step_temp(1,k) Step_Height(k+1:end)];
else
    Step_Height=[Step_Height(1:k-1) step_temp(3,k)+(-pi/2)*step_temp(1,k) step_temp(3,k)+pi/2*step_temp(1,k)];
end
end

end
% sort all the step and set teh height for the steps
xstep=sort([Step_Left Step_Right]);
xstep_height=zeros(1, 2*length(Step_Height));
xstep_height(1:2:2*length(Step_Height)-1)=Step_Height;
xstep_height(2:2:2*length(Step_Height))=Step_Height;


%plot 
figure;
plot(NumVector,Measured_Y,'o--','color',[0.8 0.8 0.8]);
label_step=ceil(len/6);% the interbal of x axis label for off periods.
set(gca,'xtick',[]);
set(gca,'xtick',NumVector(1:label_step:end));
set(gca,'xticklabel',{round(realtime_avi_rec(1:label_step:end))});  

hold on;
plot(xstep,xstep_height,'color','k','LineWidth',2);
%plot(xstep,xstep_height,'color','k','LineWidth',2);
xlabel('Time (s)');
ylabel('Frequency (Hz)');

xlim([min(NumVector) max(NumVector)]);
box on;
