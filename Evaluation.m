% Evaluation of fitting 
point=[];
counter_position=[];
counter_positionl=[];
point=[xl len];
step_dur=diff(point);
best_position=xl(2:end);
best_height=step_height;
counter_position=zeros(1,length(best_position));
counter_height=zeros(1,length(best_height));
for i = 1: length(counter_position)
counter_position(i)=best_position(i)+randi([-min(abs(step_dur)) min(abs(step_dur))],1,1);
end
for i=1 :length(counter_height)
counter_height(i)=best_height(i)+rand(1,1);   
end
counter_positionl=counter_position-1;
counter_position=[1 counter_position];
counter_positionl=[counter_positionl len];
xstep=sort([counter_position counter_positionl]);
xstep_height=zeros(1, 2*length(counter_height));
xstep_height(1:2:2*length(counter_height)-1)=counter_height;
xstep_height(2:2:2*length(counter_height))=counter_height;

ydata(1,:)=real_peak_off;
for i=1: length(counter_height)
ydata(2,counter_position(i):counter_positionl(i))=counter_height(i);
end

for i=1: length(best_height)
ydata(3,xl(i):xu(i))=best_height(i);
end

chi_best= sum((ydata(3,:)-ydata(1,:)).^2);
chi_counter= sum((ydata(2,:)-ydata(1,:)).^2);
S=chi_counter/chi_best;
