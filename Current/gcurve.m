function [ramp] = gcurve(deg,winsize,duty,delay)


slope=deg*winsize;
up=logspace(-2.8,1,slope)/10;
down=logspace(1,-2.8,slope)/10;
mid=ones(1,(round(duty*winsize)-(round(2*slope,0))));
z=zeros(1,winsize);
step=[up,mid,down];

ramp= z+[step,zeros(1,winsize-length(step))];
%hold on
%plot(ramp,'b');axis([-100 4100 -0.1 1.1])
%hold off
%end

