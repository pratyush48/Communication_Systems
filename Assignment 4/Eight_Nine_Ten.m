clear all;

close all;

[y,fs] = audioread('whatareyou2.wav');

tester = y(:,1);

freq=-fs/2:fs/(length(tester)-1):fs/2;

t = 2/fs:1/fs:9.0413;
figure(1);
plot(freq, abs(fftshift(fft(tester))));
title("Original");
% figure(2);
% plot(freq, abs(fftshift(fft(tester1))));
% title("Guassian");
% figure(3);
% plot(freq, abs(fftshift(fft(tester3))));
% title("Colored Noise");
% figure(4);
% plot(freq, abs(fftshift(fft(tester4))));
% title("Impulse Noise"); 

grid

 

fc1 = 8000;

AmpMod1 = tester'.*cos(2*pi*fc1*t);
tester1 = awgn(AmpMod1,3, 'measured');
tester2 = dsp.ColoredNoise('Color','pink','SamplesPerFrame',length(y));
cn = tester2();
tester3 = 0.2*(cn)*(rms(y)/rms(cn))+AmpMod1.';
ind = randi([0,2],1,length(y));
ind = 0.8*ind*(rms(y)/rms(ind));
tester4 = AmpMod1.'+ind.';
figure

plot(freq, abs(fftshift(fft(AmpMod1))));

title("Amplitude Modulated, fc = 8000");

 

grid();

 

fc2 = 12000;

AmpMod2 = tester'.*cos(2*pi*fc2*t);

figure();

plot(freq, abs(fftshift(fft(AmpMod2))));

title("Amplitude Modulated, fc = 12000");

 

grid

 

% Demodule each using appropriately synced carrier frequency

%-------------------------------fc1---------------------------------- 
theta = pi/2;
df = 1000;
Rcv1 = AmpMod1.*cos(2*pi*(fc1)*t);
Rcv2 = AmpMod1.*cos(2*pi*(fc1 + df)*t);
Rcv3 = AmpMod1.*cos(2*pi*(fc1)*t + theta);
Rcvawgn = tester1.*cos(2*pi*(fc1)*t);
Rcvcnoise = (tester3.').*cos(2*pi*(fc1)*t);
Rcvimpnoise = (tester4.').*cos(2*pi*(fc1)*t);
figure

plot(freq, abs(fftshift(fft(Rcv1))));

title("Proper Demodulation");

grid
figure

plot(freq, abs(fftshift(fft(Rcv2))));

title("Demodulation with df");

grid
figure

plot(freq, abs(fftshift(fft(Rcv3))));

title("Demodulation with theta");

grid
figure

plot(freq, abs(fftshift(fft(Rcvawgn))));

title("Demodulation_awgn");

grid
figure

plot(freq, abs(fftshift(fft(Rcvcnoise))));

title("Demodulation_cnoise");

grid
figure

plot(freq, abs(fftshift(fft(Rcvimpnoise))));

title("Demodulation_impnoise");

grid
rcv1 = lowpass(Rcv1,5000,fs);
rcv2 = lowpass(Rcv2,5000,fs);
rcv3 = lowpass(Rcv3,5000,fs);
rcvawgn = lowpass(Rcvawgn,5000,fs);
rcvcnoise = lowpass(Rcvcnoise,5000,fs);
rcvimpnoise = lowpass(Rcvimpnoise,5000,fs);
figure

plot(freq, abs(fftshift(fft(rcv1))));

title("After passing through LPF");

grid
figure

plot(freq, abs(fftshift(fft(rcv2))));

title("After passing through LPF df");

grid

figure

plot(freq, abs(fftshift(fft(rcv3))));

title("After passing through LPF theta");

grid
figure

plot(freq, abs(fftshift(fft(rcvawgn))));

title("After passing through LPFawgn");

grid
figure

plot(freq, abs(fftshift(fft(rcvcnoise))));

title("After passing through LPFcnoise");

grid

figure

plot(freq, abs(fftshift(fft(rcvimpnoise))));

title("After passing through LPFimpnoise");

grid
rcvp1 = audioplayer(rcv1,fs);
rcvp2 = audioplayer(rcv2,fs);
rcvp3 = audioplayer(rcv3,fs);
% play(rcvp1);
%play(rcvp2);
 play(rcvp3);
%---------------------------------------fc1----------------------------