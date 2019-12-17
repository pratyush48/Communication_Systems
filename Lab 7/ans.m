oversampling_factor = 16;
%for a pulse with amplitude one, the max frequency deviation is given
% by kf
% kf = 0.25;
kf = 4;
%increase the oversampling factor if kf (and hence frequency
%deviation, and hence bw of FM s
oversampling_factor = ceil(max(kf,1)*oversampling_factor);
ts=1/oversampling_factor;%sampling time
nsamples = ceil(1/ts);
pulse = ones(nsamples,1); %rectangular pulse
nsymbols = 10;
symbols = zeros(nsymbols,1);
%random symbol sequence
symbols = sign(rand(nsymbols,1)-0.5);
%generate digitally modulated message
nsymbols_upsampled=1+(nsymbols-1)*nsamples;
symbols_upsampled=zeros(nsymbols_upsampled,1);
symbols_upsampled(1:nsamples:nsymbols_upsampled)=symbols;
message = conv(symbols_upsampled,pulse);
%FM signal phase obtained by integrating the message
theta = 2*pi*kf*ts*cumsum(message);
cenvelope=exp(j*theta);
% phi = 2*pi*rand; %phase uniform over [0,2 pi]
% cenvelope = cenvelope*exp(j*phi);
% phi = 2*pi*rand; %phase uniform over [0,2 pi]
% df = 0.3;
% cenvelope = cenvelope.*exp(j*(2*pi*df*time+phi));
L=length(cenvelope);
time=(0:L-1)*ts;
pulse = sin(pi*time).*ones(nsamples,1);
Icomponent = real(cenvelope);
Qcomponent= imag(cenvelope);
%plot I component................................
figure(1);
plot(time,Icomponent);
figure(2);
plot(time,Qcomponent);
figure(3);
plot(time, theta/pi)
%plot a2.........................................
figure(4);
plot(time, Icomponent);
figure(5);
plot(time, Qcomponent);
%plot a3..........................................
%baseband discriminator
%differencing operation approximates derivative
Iderivative = [0;diff(Icomponent)]/ts;
Qderivative = [0;diff(Qcomponent)]/ts;
dtheta = (Icomponent.*Qderivative - Qcomponent.*Iderivative)/(Icomponent.*Icomponent + Qcomponent.*Qcomponent);
m_est = (1/(2*pi*kf)*dtheta);
figure(6);
subplot(2,1,1);
plot(time, message);
subplot(2,1,2);
plot(time, m_est);
[F, freq] = FFT(cenvelope);
figure(7);
plot(freq, abs(F));
%.......................................................................
nsymbols =1000;
symbols=zeros(nsymbols,1);
nruns=1000; 
fs=0.1;
fs_desired = 0.1;
Nmin = ceil(1/(fs_desired*ts)); %minimum length DFT for desired frequency granularity
message_length=1+(nsymbols-1)*nsamples+length(pulse)-1;
Nmin = max(message_length,Nmin);
% %for efficient computation, choose FFT size to be power of 2
Nfft = 2^(nextpow2(Nmin)); %FFT size = the next power of 2 at least as big as Nmin
psd=zeros(Nfft,1);
for runs=1:nruns,
%random symbol sequence
symbols = sign(rand(nsymbols,1)-0.5);
nsymbols_upsampled=1+(nsymbols-1)*nsamples;
symbols_upsampled=zeros(nsymbols_upsampled,1);
symbols_upsampled(1:nsamples:nsymbols_upsampled)=symbols;
% disp(pulse);
% disp("hello");
% disp(length(symbols_upsampled));
message = conv(symbols_upsampled,pulse);
%FM signal phase
theta = 2*pi*kf*ts*cumsum(message);
cenvelope=exp(j*theta);
time=(0:length(cenvelope)-1)*ts;
% %freq domain signal computed using DFT
cenvelope_freq = ts*fft(cenvelope,Nfft); %FFT of size Nfft, automatically zeropads as needed
cenvelope_freq_centered = fftshift(cenvelope_freq); %shifts DC to center of spectrum
psd=psd+abs(cenvelope_freq_centered).^2;
end
psd=psd/(nruns*nsymbols);
fs=1/(Nfft*ts); %actual frequency resolution attained
% %set of frequencies for which Fourier transform has been computed using DFT
freqs = ((1:Nfft)-1-Nfft/2)*fs;
%plot the PSD
figure(8);
disp(psd);
plot(freqs,psd); 
function [signal_freqdomain_centered,freqs] = FFT(u)
    ts = 1/16;
    fs_desired = 1/10; 

    Nmin = ceil(1/(fs_desired*ts)); 

    Nfft = 2^(nextpow2(Nmin)) ;
    signal_freqdomain = ts*fft(u,Nfft);
    signal_freqdomain_centered = fftshift(signal_freqdomain);

    fs=1/(Nfft*ts); 

    freqs = ((1:Nfft)-1-Nfft/2)*fs;
end
