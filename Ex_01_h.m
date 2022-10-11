clc
clear
close all


T = 100;
tau = 100;
zeta = 0.25;
kp = 1.1219e-05;
tmin = 0;
tmax = 1000;
delta_t = 1;
t = [tmin:delta_t:tmax];
Q = [ zeros(1,max(size([tmin:delta_t:T]))),...
ones(1,max(size([T+delta_t:delta_t:tmax]))) ];
[val,tindex] = min(abs(t-(T+delta_t)));
T2 = [ zeros(1,max(size([tmin:delta_t:T]))) ];
for ii=tindex:max(size(Q))
T2_next = inv(1/delta_t^2 + 2/delta_t + 2)*...
[kp*Q(ii) + ...
(2/delta_t^2+2/delta_t)*T2(ii-1) ...
- (2/delta_t^2)*T2(ii-2)];
T2 = [T2 T2_next];
end
plot(t,T2)