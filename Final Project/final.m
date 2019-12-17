clear all; 
close all; 
filename = 'kaaki.wav';
inter = 'whatareyou2.wav';
[err,er_fs] = audioread(inter);
[new_N,new_C] = size(inter); 
new_ts=1/er_fs;
t_er=(0:new_ts:(new_N*new_ts)- new_ts)';
fer = 10000;
%Create the message signal
% f1=50;%Modulating frequency
% msg=sin(2*pi*f1*t);
msg = err;
kf=.0628;%Modulation index
%Create the real and imaginary parts of a CW modulated carrier to be tracked.
noise =exp(1j*(2*pi*fer*t_er+2*pi*kf*cumsum(msg)));
[x,fs] = audioread(filename);
sound(x,fs);
[t,Signal,fs,N,C,new_t] = fm_transmitter(x,fs);
for i = 1:1:length(Signal)
Signal(i) = Signal(i) + noise(i);
end
Signal = awgn(Signal,25);
demod = fm_receiver(Signal,fs);

%Initilize PLL Loop 
sound(demod,fs);

function [demodulated] = fm_receiver(Signal,fs)
rng('default');
fc_arr = randi([11000 12000],1,6)
% fc_arr = r*1000;
k = size(Signal);
signal_re = buffer(Signal,ceil(k(1)/6));
demodulated = [];
    for i = 1:6
        fc = fc_arr(i);
        signal = signal_re(:,i);
        phi_hat(1)=30; 
        e(1)=0; 
        phd_output(1)=0; 
        vco(1)=0; 
        %Define Loop Filter parameters(Sets damping)
        kp=0.15; %Proportional constant 
        ki=0.1; %Integrator constant 
        %PLL implementation 
        for n=2:length(signal) 
        vco(n)=conj(exp(1j*(2*pi*n*fc/fs+phi_hat(n-1))));%Compute VCO 
        phd_output(n)=imag(signal(n)*vco(n));%Complex multiply VCO x Signal input 
        e(n)=e(n-1)+(kp+ki)*phd_output(n)-ki*phd_output(n-1);%Filter integrator 
        phi_hat(n)=phi_hat(n-1)+e(n);%Update VCO 
        end
        demodulated = [demodulated e];
    end
end

function [t,Signal,fs,N,C,new_t] = fm_transmitter(x,fs)

f=12000;%Carrier frequency 
[N,C] = size(x);
Ts=1/fs;
t=(0:Ts:(N*Ts)- Ts);
new_x = buffer(x,ceil(N/6));
t = transpose(t);
new_t = buffer(t,ceil(N/6));
% [N,C] = size(x);

% t=(0:Ts:(N*Ts)- Ts);
%Create the message signal
% msg = x';
kf=.228;%Modulation index
%Create the real and imaginary parts of a CW modulated carrier to be tracked.
rng('default');
arr = randi([11000 12000],1,6);
% arr = [1,4,5,3,2,6]*2000;
Signal = [];
for i = 1:6
    Signal= [Signal;exp(1j*(2*pi*arr(i)*new_t(:,i)+2*pi*kf*cumsum(new_x(:,i))))];%Modulated carrier
end
end