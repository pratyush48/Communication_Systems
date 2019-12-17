dt = 0.001;
to = 2;
theta = pi/4;
t = -7:dt:7;
u = signalx(t);
v = signalx1(t);
s = u + 1i*v;
sc = signalx(-t) - 1i*signalx1(-t);
figure(3);
plot(t,double(real(s)));
hold on;
plot(t,double(real(sc)));
hold off;
xlabel({'Time','(In Seconds)'});ylabel({'Amplitude','(In Units)'});
title('Real parts of s(t) and smf(t)');
legend('s(t)','smf(t)');
grid on;
figure(4);
plot(t,imag(s));
hold on;
plot(t,imag(sc));
hold off;
xlabel({'Time','(In Seconds)'});ylabel({'Amplitude','(In Units)'});
title('Imaginary parts of s(t) and smf(t)');
legend('s(t)','smf(t)');
grid on;
[s_c,t1] = contconv(double(s),double(sc),t(1),t(1),dt);
figure(5);
plot(t1,real(s_c));
xlabel({'Time','(In Seconds)'});ylabel({'Amplitude','(In Units)'});
title('Real part of convulution between s(t) and smf(t)');
legend('real(s(t)*smf(t))');
grid on;
figure(6);
plot(t1,imag(s_c));
xlabel({'Time','(In Seconds)'});ylabel({'Amplitude','(In Units)'});
title('Imaginary part of convulution between s(t) and smf(t)');
legend('imag(s(t)*smf(t))');
grid on;
figure(7);
plot(t1,abs(s_c));
xlabel({'Time','(In Seconds)'});ylabel({'Amplitude','(In Units)'});
title('Magnitude of convulution between s(t) and smf(t)');
legend('abs/mag(s(t)*smf(t))');
grid on;
s1 = (signalx(t-2) + 1i*signalx1(t-2))*exp(1i*theta);
[s1_c,t3] = contconv(double(s1),double(sc),t(1),t(1),dt);
figure(8);
plot(t3,real(s1_c));
xlabel({'Time','(In Seconds)'});ylabel({'Amplitude','(In Units)'});
title('Real part of convulution between s1(t) and smf(t)');
legend('real(s1(t)*smf(t))');
grid on;
figure(9);
plot(t3,imag(s1_c));
xlabel({'Time','(In Seconds)'});ylabel({'Amplitude','(In Units)'});
title('Imaginary part of convulution between s1(t) and smf(t)');
legend('imag(s1(t)*smf(t))');
grid on;
figure(10);
plot(t3,abs(s1_c));
xlabel({'Time','(In Seconds)'});ylabel({'Amplitude','(In Units)'});
title('Magnitude of convulution between s1(t) and smf(t)');
legend('abs/mag(s1(t)*smf(t))');
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel',xt/1)
grid on;
grid minor;
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
function [y,t] = contconv(x1,x2,s1,s2,dt)
y = conv(x1,x2)*dt;
s1_2 = s1 + (length(x1)-1)*dt;
s2_2 = s2 + (length(x2)-1)*dt;
t1 = s1+ s2;
t2 = s2_2 + s1_2;
t = t1:dt:t2;
end
