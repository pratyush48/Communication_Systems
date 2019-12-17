
ns = 4;

[m_t, time_m] = signalm(ns);

figure(1);
plot(time_m, m_t);

[u_t, time_u] = DSB(m_t, time_m, 10);
[U_f, freq] = PSD(u_t, time_u(1), 4 * ns);

figure(2);
plot(time_u, u_t);
for i = 1:2000
    [U_f_new, freq_new] = PSD(u_t, time_u(1), 4 * ns);
    U_f = add_to_avg(U_f, i, U_f_new);
end
figure(3);
plot(freq, abs(U_f) .* abs(U_f)/4);

function [U_f_avg] = add_to_avg(U_f_in,  times, U_f_new)
    
    U_f_avg = (U_f_in * times + U_f_new)/(times+1);

end

function [f, t] = DSB(m_t, time_m, fc)

f = m_t .* (cos(2 * pi * fc * time_m)); 
t = time_m;

end
function [f, freq] = PSD(m_t, time_m_start, m)

[U_f, freq, ~] = contFT(m_t, time_m_start, 1/m, 1000);
f = abs(U_f);

end

function [X, f, df] = contFT(x, tstart, dt, df_desired)
    Nmin = max(ceil(1/(df_desired * dt)), length(x));
    Nfft = 2^(nextpow2(Nmin));
    X = dt * fftshift(fft(x, Nfft));
    df = 1/(Nfft * dt);
    f = ((0: Nfft - 1) - Nfft/2) * df;
    X = X .* exp(-1i * 2 * pi * f * tstart);
end

function f = random_symbols(ns)
    A = [-1, 1];
    symbols = 1:ns;
    for i = 1:ns
        pos = randi(2);
        symbols(i) = A(pos);
    end
    f = symbols;
end

function [f, t] = signalm(ns)
    
m=4 * ns; %sampling rate as multiple of symbol rate
%discrete time representation of sine pulse
time_p = 0:1/m:4; %sampling times over duration of pulse
p = sin(pi*time_p); %samples of the pulse
%symbols to be modulated

% symbols = random_symbols(ns);
symbols = [-1;1;1;-1];
%UPSAMPLE BY m
nsymbols = length(symbols);%length of original symbol sequence
nsymbols_upsampled = 1+(nsymbols-1)*m;%length of upsampled symbol sequence
symbols_upsampled = zeros(nsymbols_upsampled,1);%
symbols_upsampled(1:m:nsymbols_upsampled)=symbols;%insert symbols with spacing M
%GENERATE MODULATED SIGNAL BY DISCRETE TIME CONVOLUTION
f=conv(symbols_upsampled,p);
%PLOT MODULATED SIGNAL
t = 0:1/m:(length(f)-1)/m; %unit of time = symbol time T

end