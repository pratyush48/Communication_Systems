[U_f, freq] = PSD(4);

% for i = 1:2000
%     [U_f_new, freq_new] = PSD(4);
%     U_f = add_to_avg(U_f, i, U_f_new);
% end
figure(2);
plot(freq, abs(U_f) .* abs(U_f)/4);

function [U_f_avg] = add_to_avg(U_f_in,  times, U_f_new)
    
    U_f_avg = (U_f_in * times + U_f_new)/(times+1);

end

function [f, freq] = PSD(ns)
m=4 * ns; %sampling rate as multiple of symbol rate
%discrete time representation of sine pulse
time_p = 0:1/m:ns; %sampling times over duration of pulse
p = sin(pi*time_p); %samples of the pulse
%symbols to be modulated

symbols = random_symbols(ns);

%UPSAMPLE BY m
nsymbols = length(symbols);%length of original symbol sequence
nsymbols_upsampled = 1+(nsymbols-1)*m;%length of upsampled symbol sequence
symbols_upsampled = zeros(nsymbols_upsampled,1);%
symbols_upsampled(1:m:nsymbols_upsampled)=symbols;%insert symbols with spacing M
%GENERATE MODULATED SIGNAL BY DISCRETE TIME CONVOLUTION
u=conv(symbols_upsampled,p);
%PLOT MODULATED SIGNAL
time_u = 0:1/m:(length(u)-1)/m; %unit of time = symbol time T
% figure(1);

% xlabel("t/T");
fc = 10;

u = u .* (cos(2 * pi * fc * time_u));
plot(time_u,u);
[U_f, freq, ~] = contFT(u, time_u(1), 1/m, 1000);
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
