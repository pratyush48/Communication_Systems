function plot_two
dt = 0.01;%sample spacing 
s1 = -2:dt:-1; %sampling times over the interval [-2,-1]
s2 = 1:dt:3; %sampling times over the interval [1,3]
x1 = 3*ones(length(s1),1); %samples for first box
x2 = 4*ones(length(s2),1); %samples for second box
[y,t] = contconv(x1,x2,s1(1),s2(1),dt);
[y1,t1] = contconv(x1,x1,s2(1),s2(1),dt);
figure(1);
subplot(2,2,[3,4]);plot(t,y);xlabel({'Time','(In Seconds)'});ylabel({'Amplitude','(In Units)'});
title('Convolution of 3I[-2,-1] and 4I[1,3]');
legend('conv(x1,x2)');
grid on;
subplot(2,2,1);plot(s1,x1);xlabel({'Time','(In Seconds)'});ylabel({'Amplitude','(In Units)'});
title('3I[-2,-1]');
grid on;
legend('x1(t)');
subplot(2,2,2);plot(s2,x2);xlabel({'Time','(In Seconds)'});ylabel({'Amplitude','(In Units)'});
title('4I[1,3]');
grid on;
legend('x2(t)');
figure(2);
subplot(2,1,2);plot(t1,y1);xlabel({'Time','(In Seconds)'});ylabel({'Amplitude','(In Units)'});
title('Convolution of 3I[-2,-1] and 3I[-2,-1]');
grid on;
legend('conv(x1,x1)');
subplot(2,1,1);plot(s1,x1);xlabel({'Time','(In Seconds)'});ylabel({'Amplitude','(In Units)'});
title('3I[-2,-1]');
grid on;
legend('x1(t)');

end
function [y,t] = contconv(x1,x2,s1,s2,dt)
y = conv(x1,x2)*dt;
s1_2 = s1 + (length(x1)-1)*dt;
s2_2 = s2 + (length(x2)-1)*dt;
t1 = s1+ s2;
t2 = s2_2 + s1_2;
t = t1:dt:t2;
end

