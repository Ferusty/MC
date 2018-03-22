function [ramp] = grwin(deg,winsize,duty,delay)



slope=deg*winsize;
up=linspace(0,1,slope);
down=linspace(1,0,slope);
mid=ones(1,(round(duty*winsize)-(round(2*slope,0))));
z=zeros(1,winsize);
step=[up,mid,down];
nil=[delay step];
ramp= z+[nil,zeros(1,winsize-length(nil))];

%hold on
%plot(ramp,'b');axis([-100 4100 -0.1 1.1])
%hold off

end

