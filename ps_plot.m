function y= ps_plot(spec_sampling,periods,peak,real_peak,pottsL2,step_modified,freq_int,Inter,real_time,label_step,minfre,maxfre,file_name,rec)

        Inter=Inter./spec_sampling;
        pottsL3 = pottsL2;
        z=10*ones(1,length(peak));
        hold on;
        timescale = linspace(1,length(peak),length(peak));
        plot3(timescale,real_peak,z,'.','color','k','MarkerSize',8);
        plot3(timescale,pottsL3,z,'g--.', 'MarkerSize', 10);
        plot3(timescale,step_modified,z,'r--.', 'MarkerSize', 10);

%         plot lines;
        yline=[min(freq_int) max(freq_int)];
        zline=[20 20];
        xline=[0 0];
hold on; 
if rec==0
for i=1:length(Inter)
    plot3(xline+Inter(i),yline,zline,'--w');
end
end
  
  

  axis tight;       
%set x labels
  hold on;
  set(gca,'xtick',[]);
  set(gca,'xtick',Inter(1:label_step:end));
  set(gca,'xticklabel',{floor(periods(1:label_step:end,1))}); 
  legend('Intensity of WT','Peak','Step','Step modified','Location','best');
  
end
