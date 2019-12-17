function three()
    
    dt = 0.01;
    s1 = -2: dt: 2;
    
%     the signal p(t)
    p_t = double(signalp(s1));
    
%     the symbols bc[n] and bs[n]
    b_c_n = [-1, 1, -1, 1, -1, 1, -1, 1, -1, 1];
    b_s_n = [1, -1, 1, -1, 1, -1, 1, -1, 1, -1];
    
%     upsampling bc[n] to have the same sampling frequency as p(t)
    b_c_n_upsampled = upsample(b_c_n, 1/dt);
    
%     convolving upsampled bc[n] with p(t) to give uc(t)
    [u_c_t, time_u_c] = contconv(b_c_n_upsampled, p_t, 0, s1, dt);
    u_c_t = u_c_t.';
    figure(1);
    plot(time_u_c, u_c_t);
    xlim([0, 13]);
    ylim([-2, 2]);
    grid on;
    grid minor;
    xlabel("Time (s)");
    ylabel("Amplitude (units)");
    title("Plot of uc(t) for 10 symbols");
    
%     upsampling bs[n] to have the same sampling frequency as p(t)
    b_s_n_upsampled = upsample(b_s_n, 1/dt);
     
    %     convolving upsampled bs[n] with p(t) to give us(t)
    [u_s_t, time_u_s] = contconv(b_s_n_upsampled, p_t, s1(1), s1(1), dt);
    u_s_t = u_s_t.';
    figure(2);
    plot(time_u_s, u_s_t);
    xlim([0, 13]);
    ylim([-2, 2]);
    grid on;
    grid minor;
    xlabel("Time (s)");
    ylabel("Amplitude (units)");
    title("Plot of us(t) for 10 symbols");
    
%     multiplying uc(t) by multiplying cos(40?t) to give up1(t) 
    u_p_1_t = double(u_c_t .* cos(time_u_c .* 40 * pi));
    figure(3);
    plot(time_u_c, u_p_1_t);
    xlim([2, 6]);
    ylim([-1.5, 1.5]);
    grid on;
    grid minor;
    xlabel("Time (s)");
    ylabel("Amplitude (units)");
    title("Plot of upconverted signal up1(t)");
    
%     forming up(t) by taking the inphase and quadrature component, viz
%     uc(t)*cos(40?t) - us(t)*sin(40?t)
    u_p_t = double(u_p_1_t - u_s_t .* sin(time_u_s .* 40 * pi));
    time_u_p = time_u_c;
    figure(4);
    plot(time_u_p, u_p_t);
    xlim([2, 6]);
    grid on;
    grid minor;
    xlabel("Time (s)");
    ylabel("Amplitude (units)");
    title("Plot of upconverted signal along with in-phase and quadrature components up(t)");
    
%     The lowpassfilter h(t)
    h_t = double(signalh(s1));
    
%     reconstructing uc(t) by convolving 2*up(t)*cos(40?t + ?) with the
%     lowpass filter where ? = 0
    [v_c_t, time_v_c] = contconv(double(2 * u_p_t .* cos(time_u_p .* 40 * pi)), h_t, time_u_p(1), s1(1), dt);
    figure(5);
    plot(time_v_c, v_c_t * dt);
    xlim([0 16]);
    grid on;
    grid minor;
    xlabel("Time (s)");
    ylabel("Amplitude (units)");
    title("Plot of reconstructed inphase component of up(t) with theta = 0");
    
%     reconstructing uc(t) by convolving 2*up(t)*cos(40?t + ?) with the
%     lowpass filter where ? = 0
    [v_s_t, time_v_s] = contconv(double(2 * u_p_t .* sin(time_u_p .* 40 * pi)), h_t, time_u_p(1), s1(1), dt);
    figure(6);
    plot(time_v_s, v_s_t * dt);
    xlim([0 16]);
    grid on;
    grid minor;
    xlabel("Time (s)");
    ylabel("Amplitude (units)");
    title("Plot of reconstructed quadrature component of up(t) with theta = 0");

%     reconstructing uc(t) by convolving 2*up(t)*cos(40?t + ?) with the
%     lowpass filter where ? = /4
    [v_c_t_1, time_v_c_1] = contconv(double(2 * u_p_t .* cos((time_u_p .* 40 * pi) + (pi/4))), h_t, time_u_p(1), s1(1), dt);
    figure(7);
    plot(time_v_c_1, v_c_t_1 * dt);
    xlim([0 16]);
    grid on;
    grid minor;
    xlabel("Time (s)");
    ylabel("Amplitude (units)");
    title("Plot of reconstructed inphase component of up(t) with theta = pi/4");

%     reconstructing uc(t) by convolving 2*up(t)*sin(40?t + ?) with the
%     lowpass filter where ? = /4
    [v_s_t_1, time_v_s_1] = contconv(double(2 * u_p_t .* sin((time_u_p .* 40 * pi) + (pi/4))), h_t, time_u_p(1), s1(1), dt);
    figure(8);
    plot(time_v_s_1, v_s_t_1 * dt);
    xlim([0 16]);
    grid on;
    grid minor;
    xlabel("Time (s)");
    ylabel("Amplitude (units)");
    title("Plot of reconstructed quadrature component of up(t) with theta = pi/4");
    
%     reconstructing given theta and vs(t)and vc(t)
    u_c_t_1 = (v_c_t_1 * cos(pi/4) + v_s_t_1 * sin(pi/4)) / 25;
    u_s_t_1 = (v_c_t_1 * sin(pi/4) - v_s_t_1 * sin(pi/4)) /25 ;
    figure(9);
    plot(time_v_s_1, u_c_t_1);
    xlim([0 16]);
    grid on;
    grid minor;
    xlabel("Time (s)");
    ylabel("Amplitude (units)");
    title("Plot of uc(t) after reconstruction from given theta and vc(t) and vs(t)");
    
    figure(10);
    plot(time_v_s_1, u_s_t_1);
    xlim([0 16]);
    grid on;
    grid minor;
    xlabel("Time (s)");
    ylabel("Amplitude (units)");
    title("Plot of us(t) after reconstruction from given theta and vc(t) and vs(t)");
    
    
end

%the signal p(t)
function f = signalp(time_steps)
    %Creating a symbolic function using piecewise
    syms x;
    out = piecewise(0 <= x <= 1, 1, 0);    
    % Substituting the time steps to the symbolic function
    f = subs(out, x, time_steps); 
end

function f = signalh(time_steps)
    %Creating a symbolic function using piecewise
    syms x;
    out = piecewise(0 <= x <= 0.25, 1, 0);    
    % Substituting the time steps to the symbolic function
    f = subs(out, x, time_steps); 
end

%calculates convolution of two given signals
function [y, t] = contconv(x1, x2, t1, t2, dt)
    y = conv(x1, x2);
    t = 0: dt: ((length(y) -1)*dt);
end

%upsamples the given function 
function f = upsample(x1, m)
    nsymbols = length(x1);
    nsymbols_upsampled = 1 + (nsymbols - 1) * m;
    symbols_upsampled = zeros(nsymbols_upsampled, 1);
    symbols_upsampled(1: m: nsymbols_upsampled) = x1;
    
    f = symbols_upsampled;

end