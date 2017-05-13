function stepvar=stepfit_test(x,y,len)
%x is the only the peak number, here don't set the real time.
%y is the measured data need to be fitted.
Nw=length(x);% the step window (number of x in the plateau)

%set three parameters for the fittype: actan(kx);
k1=1e3/Nw;
k2=k1*10;
k3=k2*10;

%Initial point, use the number of elements in, and size of, x0 to 
%determine the number and size of variables that fun accepts.
x0 = [(y(end)-y(1))/2,mean(x),mean(y)];
%Lower bounds of initial value;
lb=[-10,0, 0];
%Upper bounds of initial value;
lu=[10,len, max(y)];
%Set the parameters for options in the fitting:
%TolFun: Termination tolerance on the function value
%MaxFunEvals: Maximum number of function evaluations allowed.
%MaxIter: Maximum number of iterations allowed.
options = optimset('TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',1e6,'MaxIter',1e6);
% Create the first fitting function
fun = @(xx,xdata)xx(1)*atan(k1*(xdata-xx(2)))+xx(3);
% first fitting
xx1 = lsqcurvefit(fun,x0,x,y,lb,lu,options);
%plot(x,y,'ro',x,fun(xx1,x),'r-')
% Create the second fitting function
fun = @(xx,xdata)xx(1)*atan(k2*(xdata-xx(2)))+xx(3);
% second fitting
xx2 = lsqcurvefit(fun,xx1,x,y,lb,lu,options);
%plot(x,y,'go',x,fun(xx2,x),'k-')
% Create the third fitting function
fun = @(xx,xdata)xx(1)*atan(k3*(xdata-xx(2)))+xx(3);
% third fitting
stepvar = lsqcurvefit(fun,xx2,x,y,lb,lu,options);
 %plot(x,y,'o--','color',[0.8 0.8 0.8]);
%hold on;
 %plot(x,fun(stepvar,x),'k-','LineWidth',2);

end
