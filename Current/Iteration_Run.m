clear;

format longE;
trials=100;
pNMSE=[];

for i=1:trials

%% Parameters for input signal and multicoset sampler
W=1000; % W/2 Hz is bandlimit of multiband signal class (W=Nyquist rate (Hz))
B=10; % maximum bandwidth of each occupied frequency band (Hz)
K=5; % number of positive disjoint frequency bands present in input signal x
T=20; % duration of observation interval in sec
q=20; % q (uniform) samples collected every L/W seconds (each channel samples W/L Hz)
L=100; % L odd doesn't work
nyqover=1; % controls the length of the simulated multiband signal x (conceptually, x is a vector 
% formed from sampling a continuous-time multiband signal at a rate that is nyqover times more
% than this Nyquist rate. nyqover is a positive integer.) 
SNR=100; % signal-to-noise ratio (in dB)

%% Parameter check
if q>L
 error('q must be less than or equal to L to ensure sub-Nyquist sampling');
end

%% Simulated signal generation, (sub)sampling, support recovery, signal reconstruction
[x,centerfreq,t]=bandsparse(W,B,K,L,nyqover,'FH2H',T);
input_sigP=norm(x)^2/length(x); % compute signal power
input_noiseP=(1/(10^(SNR/10)))*input_sigP; % determine required noise power (variance) wrt input signal power
randn('state',8556);
%w=sqrt(input_noiseP)*(randn(s,size(x))+1i*randn(s,size(x))); % generate white Gaussian noise with required variance 
w=sqrt(input_noiseP)*(randn(size(x))+1i*randn(size(x))); % generate white Gaussian noise with required variance 
xw=x+w; % add noise

[yw,cw]=mc_sampling(xw,L,q,nyqover); % multi-coset sampling
[eigspectw,centerfreqw_hat,xw_hat,Rw,Phiw_Omega]=mc_recovery(yw,L,W,cw,t,nyqover,'LS'); % recover support and reconstruct signal

% figure % graphically display sampling pattern and Phi matrix
% subplot(221); hist(cw,L); 
% set(gca,'ytick',[],'DataAspectRatio',[L 3 1]); title('Sampling pattern');
% subplot(222); imagesc(abs(Phiw_Omega)); 
% axis equal tight; title('\Phi matrix (magnitude)');
% subplot(223); imagesc(abs(Rw)); 
% axis equal tight; title('R matrix (magnitude)');
% saveas(gcf,'Matrix.jpg')
% 
% 
% figure % plot fft of original signal x and reconstructed signal x_hat
% f= (-length(x)/2:length(x)/2-1)./length(x) * (W*nyqover);
% h0=subplot(311); plot(f,abs(fftshift(fft(x))),'Marker','none'); title('Fourier Spectrum (original signal)'); 
% xlabel('Hz'); ylabel('Magnitude'); xlim([-W/2 W/2]); 
% h1=subplot(312); plot(f,abs(fftshift(fft(xw))),'Marker','none'); title('Fourier Spectrum (noisy signal)'); 
% xlabel('Hz'); ylabel('Magnitude')
% h2=subplot(313); plot(f,abs(fftshift(fft(xw_hat))),'Marker','none'); title('Fourier Spectrum (reconstructed signal)');
% xlabel('Hz'); ylabel('Magnitude')
% linkaxes([h0 h1 h2])
% saveas(gcf,'Fourier.jpg')

% figure % display support recovery related quantities
% subplot(211); plot(abs(eigspectw),'.-'); title('Modified MUSIC Spectrum'); 
% xlabel('Column index of \Phi matrix'); axis([1 L 0 max(abs(eigspectw))])
% subplot(212); bar(-W/2+W/L:W/L:W/2,histc(centerfreq,-W/2+W/(2*L):W/L:W/2),'EdgeColor','w','FaceColor','r'); 
% hold on; bar(-W/2+W/L:W/L:W/2,histc(centerfreqw_hat,-W/2+W/(2*L):W/L:W/2),'EdgeColor','w');  
% set(gca,'ytick',[]); xlabel('Hz'); axis([-W/2 W/2 0 1]); title('Spectral slices identified as active') 
% saveas(gcf,'MUSIC Spectrum.jpg')

% figure % plot original multiband signal x and reconstructed signal x_hat
% h3=subplot(411); plot(t,real(x),'Marker','none');
% xlabel('Seconds'); title('Original signal (real part)');
% h4=subplot(412); plot(t,real(xw),'Marker','none');
% xlabel('Seconds'); title('Noisy signal (real part)');
% h5=subplot(413); plot(t,real(xw_hat),'r','Marker','none'); 
% xlabel('Seconds'); title('Reconstructed signal (real part)');
% h6=subplot(414); stem(t,real(xw-xw_hat),'Marker','none'); 
% xlabel('Seconds'); title('Difference Signal (noisy-reconstructed)');
% linkaxes([h3 h4 h5 h6])
% saveas(gcf,'Reconstruction.jpg')

% figure %plot specgram
% specgram(x);
% saveas(gcf,'Spectrogram.jpg')
% title('Spectrogram');

%discontinuity causing smearing. significant energy in frequency
%significant energy in frequency, when you expect maybe one or two spectral
%slices to be active there is a spillage.

% fprintf('Nyquist rate %3.3f Hz\n',W);
% fprintf('System sampling rate %3.3f Hz\n',(q*W)/L);
% fprintf('Subsampling factor %3.3f\n',q/L);

NMSE=calNMSE(x,xw_hat);

% fprintf('Normalized Means Square Error %2.3e \n', NMSE);
pNMSE(end+1)=NMSE;
end
stem(pNMSE);
avgNMSE=(sum(pNMSE)/length(pNMSE))

