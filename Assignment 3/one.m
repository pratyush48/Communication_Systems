function one
N = 100;
dt = 0.1;
t = -N:dt:N;
% I = signalx(t);
bc = randi([0,1],1,N);
uc = modsig(N,t,bc);
% display(I);
plot(t,uc);
end
function m = modsig(N,t,arr)
m1 = zeros(N);
for i = 1:N
    if(arr(i) == 0)
        arr(i) = -1;
    end
    m1(i) = arr(i)*signalx(t-i);
end
m = m1;
end
function u = signalx(t)
syms x;
y = piecewise(0 <= x <= 1, 1, 0);
u = subs(y,x,t);
end