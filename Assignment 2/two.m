function two
q_2();
q_3();
end
function q_2
dt = (1/16);
t = -8:dt:8;
s = 3*sinc(2*t - 3);
[Y,f,df] = contFT(s,t(1),dt,10^(-3));
figure(11);
plot(f,abs(Y));
xlabel({'Frequency','(In MHz)'});ylabel({'Amplitude','(In Units)'});
title('Magnitude of S(f)');
legend('abs/mag(S(f)))');
grid on;
figure(12);
plot(f,angle(Y));
xlabel({'Frequency','(In MHz)'});ylabel({'Amplitude','(In Units)'});
title('Phase of S(f)');
legend('angle(S(f))');
grid on;
grid minor;
end
function q_3
dt = 0.1;
to = 2;
theta = pi/4;
t = -7:dt:7;
u = signalx(t);
v = signalx1(t);
s = u + 1i*v;
sc = signalx(-t) - 1i*signalx1(-t);
[s_c,t1] = contconv(double(s),double(sc),t(1),t(1),dt);
[S,f,df] = contFT(double(s),t(1),dt,10^(-3));
[S_C,F,DF] = contFT(double(s_c),t1(1),dt,10^(-3));
figure(13);
plot(f,abs(S));
xlabel({'Frequency','(In KHz)'});ylabel({'Amplitude','(In Units)'});
title('Magnitude of |S(f)|');
legend('abs/mag(S(f))');
grid on;
figure(14);
plot(F,abs(S_C));
xlabel({'Frequency','(In KHz)'});ylabel({'Amplitude','(In Units)'});
title('Magnitude of fourier transform of convolution between s(t)  and smf(t)');
legend('abs/mag(fft(s(t)*smf(t)))');
grid on;
figure(15);
plot(F,angle(S_C));
xlabel({'Frequency','(In KHz)'});ylabel({'Amplitude','(In Units)'});
title('Phase of fourier transform of convolution between s(t) and smf(t)');
legend('angle(fft(s(t)*smf(t)))');
grid on;
end
function [y,t] = contconv(x1,x2,s1,s2,dt)
y = conv(x1,x2)*dt;
s1_2 = s1 + (length(x1)-1)*dt;
s2_2 = s2 + (length(x2)-1)*dt;
t1 = s1+ s2;
t2 = s2_2 + s1_2;
t = t1:dt:t2;
end
function u = signalx(t)
syms x;
y = piecewise(1 <= x <= 2, 2, 2 <= x <= 3, -1, 3 <= x <= 4, -3, 0);
u = subs(y,x,t);
end
function v = signalx1(t)
syms x;
y = piecewise(-1 <= x <= 0, 1, 0 <= x <= 1, 3, 1 <= x <= 2, 1, 0);
v = subs(y,x,t);
end
function [X,f,df] = contFT(x,tstart,dt,df_desired) 
%Use Matlab DFT for approximate computation of continuous time Fourier transform
%INPUTS 
%x = vector of time domain samples, assumed uniformly spaced %tstart= time at which first sample is taken
%dt = spacing between samples  
%df_desired = desired frequency resolution 
%OUTPUTS 
% X=vector of samples of Fourier transform 
%f=corresponding vector of frequencies at which samples are obtained
%df=freq resolution attained (redundant--already available from %difference of consecutive entries of f
%%%%%%%%%
%minimum FFT size determined by desired freq res or length of x 
Nmin=max(ceil(1/(df_desired*dt)),length(x));  
%choose FFT size to be the next power of 2
Nfft = 2^(nextpow2(Nmin));
%compute Fourier transform, centering around DC
X=dt*fftshift(fft(x,Nfft));
%achieved frequency resolution
df=1/(Nfft*dt); 
%range of frequencies covered 
f = ((0:Nfft-1)-Nfft/2)*df;
%same as f=-1/(2*dt):df:1/(2*dt) - df %phase shift associated with start time
X=X.*exp(-1i*2*pi*f*tstart);
end
