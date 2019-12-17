dt = 0.001;
t = -5:dt:5;
u = signalx(t);
umf = signalx(-t);
figure(1);
plot(t,u);
hold on;
plot(t,umf);
hold off;
xlabel({'Time','(In Seconds)'});ylabel({'Amplitude','(In Units)'});
title('Plot of u(t) and umf(t)');
legend('u(t)','umf(t)');
grid on;
[y,t1] = contconv(double(u),double(umf),t(1),t(1),dt);
figure(2);
plot(t1,y);
xlabel({'Time','(In Seconds)'});ylabel({'Amplitude','(In Units)'});
title('Convoluton of u(t) and umf(t)');
legend('u(t)*umf(t)');
grid on;
function u = signalx(t)
syms x;
y = piecewise(1 <= x <= 2, 2, 2 <= x <= 3, -1, 3 <= x <= 4, -3, 0);
u = subs(y,x,t);
end
function [y,t] = contconv(x1,x2,s1,s2,dt)
y = conv(x1,x2)*dt;
s1_2 = s1 + (length(x1)-1)*dt;
s2_2 = s2 + (length(x2)-1)*dt;
t1 = s1+ s2;
t2 = s2_2 + s1_2;
t = t1:dt:t2;
end