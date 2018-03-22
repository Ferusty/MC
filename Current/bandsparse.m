function [x,centerfreq,t]=bandsparse(W,B,K,L,nyqover,type,T,CF,CF2)

% Usage: [x,centerfreq,t]=bandsparse(W,B,K,L,nyqover,type,T)
% W: Nyquist frequency (bandlimit of the multiband signal= W/2)
% B: bandwidth of each occupied frequency band (Hz)
% K: number of occupied bands (positive integer)
% L: parameter that determines the period of the nonuniform sampling (every L/W seconds q samples 
%    are collected); positive integer > q
% nyqover: controls the length of the simulated multiband signal x (conceptually, x is a vector 
%          formed from sampling a continuous-time multiband signal at a rate that is nyqover times more
%          than this Nyquist rate. nyqover is a positive integer.) 
% type: specifies the type of multiband signal (options are 'modsinc' (modulated sinc waveforms), 'filterwhitenoise' (white
%       noise filtered through a notch filter), and 'chirp' (multiband chirp signal)
% T: observation interval [0,T) sec
% x: simulated multiband signal
% centerfreq: center frequencies of occupied bands
% t: vector of discrete time instances over which band sparse signal is defined


switch type

case 'modsinc'
Tend=floor(ceil(T/(1/(W*nyqover)))/L)*L/(W*nyqover);
t=(0:1/(W*nyqover):Tend-1/(W*nyqover)); % time instances over which simulated signal is defined

centerfreq=rand(1,K)*W-W/2; % randomly generate center frequencies of occupied bands
%centerfreq=[-152 20 350]; % use to manually set centerfreq; comment if randomly generating delays
%delays=rand(1,K)*T; % uncomment to randomly generate time delays of sinc pulses
delays=[2.5 8.0 15.0]; % use to manually set time offsets of sinc pulses; comment if randomly generating delays
tt=(-Tend:1/(W*nyqover):Tend-1/(W*nyqover));
x=zeros(1,length(t)); % initialize multiband signal vector

for n=1:K % generate K time-domain sinc pulses and modulate each of them to their respective center frequencies
 sincpul=sqrt(B)*sinc(B*(tt-delays(n))); % time domain sinc pulses
 x=x+sincpul(floor(Tend*W*nyqover)+1:end).*exp(1i*2*pi*centerfreq(n)*t); % modulate pulses 
end


case 'filterwhitenoise' % occupied bands are filtered white Gaussian noise 
Tend=floor(ceil(T/(1/(W*nyqover)))/L)*L/(W*nyqover);
t=(0:1/(W*nyqover):Tend-1/(W*nyqover)); % time instances over which simulated signal is defined

%centerfreq=rand(1,K)*W/2; % uncomment to randomly generate center frequencies of occupied bands
centerfreq=[20 150 350]; % use to manually set centerfreq; comment if randomly generating center frequencies
x=zeros(1,length(t)); % initialize multiband signal vector
for n=1:K
 f=[(centerfreq(n)-(B/4))/(W/2) (centerfreq(n)+(B/4))/(W/2)];
 hb=fir1(200,f,'bandpass'); % window based FIR design
 %s=RandStream.getDefaultStream; reset(s); % uncomment for repeatable results
 x=x+filter(hb,1,randn(1,length(x))); % filter white Gaussian noise
end


case 'chirp' % creates linear chirp within each occupied band
Tend=floor(ceil(T/(1/(W*nyqover)))/L)*L/(W*nyqover);
t=(0:1/(W*nyqover):Tend-1/(W*nyqover)); % time instances over which simulated signal is defined

%centerfreq=rand(1,K)*W/2; % randomly generate center frequencies of occupied bands
centerfreq=[250 350 100]; % uncomment to manually set centerfreq 
x=zeros(1,length(t)); % initialize multiband signal vector
delays=[2 8 15]; % set chirp offsets manually
durations=[2 1 2]; % set chirp durations manually
for n=1:K
 chirpsig=zeros(1,length(t));
 chirppul=chirp(t(1:durations(n)*W*nyqover),centerfreq(n),t(durations(n)*W*nyqover),centerfreq(n)+B); % generate real chirp pulse
 chirpsig(delays(n)*W*nyqover:(delays(n)+durations(n))*W*nyqover-1)=chirppul;
 x=x+chirpsig;
end


case 'FH1H' % 
Tend=floor(ceil(T/(1/(W*nyqover)))/L)*L/(W*nyqover);
t=(0:1/(W*nyqover):Tend-1/(W*nyqover)); % time instances over which simulated signal is defined

x=[]; % initialize multiband signal vector
signal1=[];


%%Frequency Generation%% Random/Manual
centerfreq= randi([0 500], 20 ,1); % Randomly Generate Frequencies
%centerfreq(:,1:20)=400; % Manually Generate Frequencies
%centerfreq = [100 200 400 600 300 100 200 100 100 300 100 200 100 100 300 100 200 100 100 300]


%%Setting the duty frequencies for individual hops
%Manual
%duty = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1];
%Fixed
dutyn=0.7;
duty(:,1:20) = dutyn;
dutyl=round((1-dutyn),2);

durations(:,1:20)=1;
fs=10000;
delay1=0;
%delay1=zeros(1,randi(dutyl*length(t)/K,1,1));
for n=1:K
    fhsssig=zeros(1,length(t));
    %t1=1/10000:1/10000:T/K*durations(n)*W*nyqover/10000;
    t1=1/fs:1/fs:(T/K*durations(n)*W*nyqover/fs);
    fhsssig=zeros(1,length(t));
    %fhsspul = sin(2*pi*10*centerfreq(n)*t1);
    fhsspul=exp(j*2*pi*10*centerfreq(n)*t1);
    %cc = 0.5*square(2*pi*0.5*t1,length(t1)/100*duty)+0.5;
    %cc = 0.5*square(2*pi*1*0.5*t1,(length(pi*t1)/200)*duty(n))+0.5;
    win1=grwin(0.05,length(pi*t1),duty(n),delay1);
    signal1 = fhsspul.*win1;
    x=[x signal1];
    
    fprintf('Frequency  %3.2f Hz\n',centerfreq(n));
end

while (length(x)~=20000)
   if length(x)<20000
    x(end+1)=1;
elseif length(x)>20000
    x(end)=[];
   end   
end


case 'FH2H' % 
Tend=floor(ceil(T/(1/(W*nyqover)))/L)*L/(W*nyqover);
t=(0:1/(W*nyqover):Tend-1/(W*nyqover)); % time instances over which simulated signal is defined

x=[]; % initialize multiband signal vector
signal1=[];
signal2=[];


%%Frequency Generation%% Random/Manual
 centerfreq= randi([-500 500], 20 ,1); % Randomly Generate Frequencies
 centerfreq2= randi([-500 500], 20 ,1); % Randomly Generate Frequencies
%  centerfreq(:,1:20)=100; % Manually Generate Frequencies
%  centerfreq2(:,1:20)=150; % Manually Generate Frequencies
%centerfreq = [100 200 100 100 300 100 200 100 100 300 100 200 100 100 300 100 200 100 100 300]


%%Setting the duty frequencies for individual hops
%Manual
%duty = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1];
%Fixed
dutyn=0.8;
duty(:,1:20) = dutyn;
dutyl=round((1-dutyn),2);

durations(:,1:20)=1;
fs=10000;

%%Setting random delay to both hoppers%%
delay1=zeros(1,randi(dutyl*length(t)/K,1,1));
delay2=zeros(1,randi(dutyl*length(t)/K,1,1));
%%Synchronize timer by setting fixed delay%%
delay1=0;
delay2=0;
%%Synchronize timer by setting delay to zero%%
%delay1=0;
%delay2=0;
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
     fprintf('Frequency Hopper 1  %3.2f Hz\n',centerfreq(n));
     fprintf('Frequency Hopper 2 %3.2f Hz\n',centerfreq2(n));
end

while (length(x)~=20000)
   if length(x)<20000
    x(end+1)=1;
elseif length(x)>20000
    x(end)=[];
   end   
end





end % switch