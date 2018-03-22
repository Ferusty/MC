W=1000; % W/2 Hz is bandlimit of multiband signal class (W=Nyquist rate (Hz))
B=10; % maximum bandwidth of each occupied frequency band (Hz)
K=5; % number of positive disjoint frequency bands present in input signal x
T=20; % duration of observation interval in sec
q=25; % q (uniform) samples collected every L/W seconds (each channel samples W/L Hz)
L=100; % L odd doesn't work
nyqover=1; % controls the length of the simulated multiband signal x (conceptually, x is a vector 
% formed from sampling a continuous-time multiband signal at a rate that is nyqover times more
% than this Nyquist rate. nyqover is a positive integer.) 
SNR=100; % signal-to-noise ratio (in dB)
CF= randi([-500 0], 20 ,1); % Randomly Generate Frequencies
CF2= randi([0 500], 20 ,1);

Tend=floor(ceil(T/(1/(W*nyqover)))/L)*L/(W*nyqover);
t=(0:1/(W*nyqover):Tend-1/(W*nyqover)); % time instances over which simulated signal is defined

x=[]; % initialize multiband signal vector
signal1=[];
signal2=[];

centerfreq=CF;
centerfreq2=CF2;
%%Frequency Generation%% Random/Manual
% centerfreq= randi([-500 0], 20 ,1); % Randomly Generate Frequencies
% centerfreq2= randi([0 500], 20 ,1); % Randomly Generate Frequencies
%  centerfreq(:,1:20)=500; % Manually Generate Frequencies
%  centerfreq2(:,1:20)=450; % Manually Generate Frequencies
%centerfreq = [100 200 100 100 300 100 200 100 100 300 100 200 100 100 300 100 200 100 100 300]


%%Setting the duty frequencies for individual hops
%Manual
%duty = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1];
%Fixed
duty(:,1:20) = 0.7;

durations(:,1:20)=1;
fs=10000;
delay1=zeros(1,randi(0.3*length(t)/K,1,1));
delay2=zeros(1,randi(0.3*length(t)/K,1,1));
for n=1:K
    fhsssig=zeros(1,length(t));
    %t1=1/10000:1/10000:T/K*durations(n)*W*nyqover/10000;
    t1=1/fs:1/fs:(T/K*durations(n)*W*nyqover/fs);
    fhsssig=zeros(1,length(t));
    %fhsspul = sin(2*pi*10*centerfreq(n)*t1);
    %fhsspul = cos(2*pi*10*centerfreq(n)*t1);
    fhsspul=exp(j*2*pi*10*centerfreq(n)*t1);
    %fhsspul2 =sin(2*pi*10*centerfreq2(n)*t1);
    %fhsspul2 =cos(2*pi*10*centerfreq2(n)*t1);
    fhsspul2=exp(j*2*pi*10*centerfreq2(n)*t1);
    %cc = 0.5*square(2*pi*0.5*t1,length(t1)/100*duty)+0.5; %Square Wave
    %cc = 0.5*square(2*pi*1*0.5*t1,(length(pi*t1)/200)*duty(n))+0.5; %SW
   
    win1=grwin(0.05,length(pi*t1),duty(n),delay1);
    win2=grwin(0.05,length(pi*t1),duty(n),delay2);
    fhsspulm1=fhsspul.*win1;
    fhsspulm2=fhsspul2.*win2;
    
    cz = fhsspulm1+fhsspulm2;
    signal1=[signal1 fhsspulm1];
    signal2=[signal2 fhsspulm2];
    x=[signal1+signal2];
%     fprintf('Frequency Hopper 1  %3.2f Hz\n',centerfreq(n));
%     fprintf('Frequency Hopper 2 %3.2f Hz\n',centerfreq2(n));
end

while (length(x)~=20000)
   if length(x)<20000
    x(end+1)=1;
elseif length(x)>20000
    x(end)=[];
   end   
end



subplot(3,1,1)
plot(t,signal1);
subplot(3,1,2)
plot(t,signal2);
subplot(3,1,3)
plot(t,x);