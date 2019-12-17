function plot_one 

%Time vector
t = -6:0.001:6;
%Output vector
yb = signalx(t);
yc = signalx(t-3);
yd = signalx(3-t);
ye = signalx(2*t);
%Plot for q1
figure(1);
plot(t,yb);
%set(gca,'XLim',[-7 7])
xlabel({'Time','(In Seconds)'});
ylabel({'Amplitude','(In Units)'});
grid on;
title('Graph of x(t) for -6 <= t <= 6');
figure(2);
plot(t,yc);
%set(gca,'XLim',[-7 7])
xlabel({'Time','(In Seconds)'});
ylabel({'Amplitude','(In Units)'});
grid on;
title('Graph for x(t - 3)');
figure(3);
plot(t,yd);
%set(gca,'XLim',[-7 7])
xlabel({'Time','(In Seconds)'});
ylabel({'Amplitude','(In Units)'});
grid on;
title('Graph for x(3-t)');
figure(4);
plot(t,ye);
%set(gca,'XLim',[-7 7])
xlabel({'Time','(In Seconds)'});
ylabel({'Amplitude','(In Units)'});
grid on;
title('Graph for  x(2t)');
end


function s = signalx(t)
s = arrayfun(@arr_signal,t);
end

function x = arr_signal(t)
if(t>=-3 && t<=-1)
    x = 2*exp(t+2);
elseif(t>=-1 && t<= 4)
    x = 2*exp(-t)*cos(2*pi*t);
else
    x = 0;
end
end