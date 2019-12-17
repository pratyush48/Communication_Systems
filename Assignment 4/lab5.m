u = dlmread('data.mat');
dt = 4/601;
t = 0:dt:4;
figure;
plot(t,u);
% m = u.*cos(2*pi*100*t);
% [Y,f,~] = contFT(m,t(1),dt,10^-3);
% figure;
% plot(f,Y);
Vd(1) = 0;

for i = 2:length(u)

    if u(i) > Vd(i-1)              % diode on (charging)

        Vd(i) = u(i);

    else                                % diode off (discharging)

        Vd(i) = Vd(i-1) - (dt/0.03)*Vd(i-1);

    end

end
% [Y,f,~] = contFT(Vd,t(1),dt,10^-3);
% figure;
% plot(f,Y);
figure;
plot(t,Vd);
u(u < 0) = 0;
% t = t(u >= 0);
h = exp(-t/0.03); % rc = 0.03
[Y,ty] = contconv(u,h,0,0,4/601);
Y1 = Y - mean(Y);
lp = lowpass(Y1,10,150);
lp1 = lowpass(Vd,10,150);
figure;
plot(ty,Y1);
xlim([0 4]);
figure;
plot(ty,lp);
xlim([0 4]);
figure;
plot(t,lp1);
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

function [y,t] = contconv(x1,x2,s1,s2,dt)
y = conv(x1,x2)*dt;
s1_2 = s1 + (length(x1)-1)*dt;
s2_2 = s2 + (length(x2)-1)*dt;
t1 = s1+ s2;
t2 = s2_2 + s1_2;
t = t1:dt:t2;
end