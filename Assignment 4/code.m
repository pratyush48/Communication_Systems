% function code
% theta = pi/4;
m =  16;
dt = 1/m;
t = -2:dt:2;
time_p = 0:1/m:1; 
p = sin(pi*time_p);
A =  [-1, 1];
symbols = 1:4;
for i = 1:4
    pos = randi(2);
    symbols(i) = A(pos);
end
 bc = [-1;1;1;-1];
[uc,tc] = upsampler(p,symbols,m,4);
figure;
plot(tc,uc);
xlabel({'Time','(In milliseconds)'});ylabel({'Amplitude','(In Units)'});
title('Message Signal');
legend('m(t)');
grid on;
[Y,f] = FFT(uc);
figure;
plot(f,((abs(Y).*abs(Y))/4));
xlabel({'Frequency','(In Hz)'});ylabel({'Amplitude','(In Units)'});
title('PSD of m(t)');
legend('S_m(f)');
Y1 = Avg(p,dt,Y);
figure;
plot(f,((abs(Y1).*abs(Y1))/4));
xlabel({'Frequency','(In Hz)'});ylabel({'Amplitude','(In Units)'});
title('Avg PSD of m(t)');
legend('S_m(f)');
grid on;
fc = 10;
m1 = 4*fc;
time_p1 = 0:1/m1:1; 
p1 = sin(pi*time_p1);
[uc1,tc1] = upsampler(p1,symbols,m1,4);
u = uc1 .* (cos(2 * pi * fc * tc1));
figure;
plot(tc1,u);
xlabel({'Time','(In milliseconds)'});ylabel({'Amplitude','(In Units)'});
title('DSB Signal');
legend('u(t)');
grid on;
[Y2,f2] = FFT(u);
Y3 = DSBAvg(p1,dt,Y2);
figure;
plot(f,((abs(Y3).*abs(Y3))/4));
xlabel({'Frequency','(In Hz)'});ylabel({'Amplitude','(In Units)'});
title('Avg PSD of u(t)');
legend('S_u(f)');
Ac = 1;
u1 = (Ac+uc1).*(cos(2 * pi * fc * tc1));
figure;
plot(tc1,u1);
xlabel({'Time','(In milliseconds)'});ylabel({'Amplitude','(In Units)'});
title('AM Signal');
legend('u(t)');
grid on;
[Y4,f4] = FFT(u1);
Y5 = AMAvg(p1,dt,Y4);
figure;
plot(f,((abs(Y5).*abs(Y5))/4));
xlabel({'Frequency','(In Hz)'});ylabel({'Amplitude','(In Units)'});
title('Avg PSD of u(t)');
legend('S_u(f)');
function Y2 = Avg(p,dt,Y)
m= 16;
for i = 1:50
    A =  [-1, 1];
    symbols = 1:4;
    for j = 1:4
        pos = randi(2);
        symbols(j) = A(pos);
    end
    [uc,tc] = upsampler(p,symbols,m,4);
    [Y1,f1] = FFT(uc);
    Y = ((Y*i)+Y1)/(i+1);
end
Y2 = Y;
end
function Y2 = DSBAvg(p,dt,Y)
fc = 10;
m= 4*fc;
symbols1 = 1:50;
for i = 1:50
    A =  [-1, 1];
    symbols = 1:4;
    for j = 1:4
        pos = randi(2);
        symbols(j) = A(pos);
    end
    [uc,tc] = upsampler(p,symbols,m,4);
    u = uc .* (cos(2 * pi * fc * tc));
    [Y1,f1] = FFT(u);
    Y = ((Y*i)+Y1)/(i+1);
end
Y2 = Y;
end
function Y2 = AMAvg(p,dt,Y)
fc = 10;
m= 4*fc;
Ac = 1;
for i = 1:50
    A =  [-1, 1];
    symbols = 1:4;
    for j = 1:4
        pos = randi(2);
        symbols(j) = A(pos);
    end
    [uc,tc] = upsampler(p,symbols,m,4);
    u1 = (Ac+uc).*(cos(2 * pi * fc * tc));
    [Y1,f1] = FFT(u1);
    Y = ((Y*i)+Y1)/(i+1);
end
Y2 = Y;
end
function [s,t] = upsampler(p,bc,m,N)
nsymbols = length(bc);
nsymbols_upsampled = 1 + (nsymbols - 1) * m;
symbols_upsampled = zeros(nsymbols_upsampled, 1);
symbols_upsampled(1:m:nsymbols_upsampled) = bc;
s = conv(symbols_upsampled,p);
t = 0:1/m:(length(s) - 1)/m;
s = s.';
end
function [signal_freqdomain_centered,freqs] = FFT(u)
    ts = 1/16;
    fs_desired = 1/1000; 

    Nmin = ceil(1/(fs_desired*ts)); 

    Nfft = 2^(nextpow2(Nmin)) ;
    signal_freqdomain = ts*fft(u,Nfft);
    signal_freqdomain_centered = fftshift(signal_freqdomain);

    fs=1/(Nfft*ts); 

    freqs = ((1:Nfft)-1-Nfft/2)*fs;
end
