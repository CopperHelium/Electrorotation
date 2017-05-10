function stepvar=stepfit(x,y)
%x is the only the peak number, here don't set the real time.
%y is the measured data need to be fitted.
Nw=length(x);% the step window (number of x in the plateau)

%set three parameters for the fittype: actan function.
k1=1e6/Nw;
k2=k1*10;
k3=k2*1;

%Initial point, use the number of elements in, and size of, x0 to 
%determine the number and size of variables that fun accepts.
x0 = [(y(end)-y(1))/2,mean(x),mean(y)];
%Lower bounds
lb=[-10,-1e3, -10];
%Upper bounds
lu=[10,1e5, 20];
%Set the parameters for options in the fitting:
%TolFun: Termination tolerance on the function value
%MaxFunEvals: Maximum number of function evaluations allowed.
%MaxIter: Maximum number of iterations allowed.
options = optimset('TolFun',1e-7,'TolX',1e-7,'MaxFunEvals',1e4,'MaxIter',1e4);
% Create the first fitting function
fun = @(xx,xdata)xx(1)*atan(k1*(xdata-xx(2)))+xx(3);
% first fitting
xx1 = lsqcurvefit(fun,x0,x,y,lb,lu,options);
% Create the second fitting function
fun = @(xx,xdata)xx(1)*atan(k2*(xdata-xx(2)))+xx(3);
% second fitting
xx2 = lsqcurvefit(fun,xx1,x,y,lb,lu,options);
% Create the third fitting function
fun = @(xx,xdata)xx(1)*atan(k3*(xdata-xx(2)))+xx(3);
% third fitting
stepvar = lsqcurvefit(fun,xx2,x,y,lb,lu,options);
end
