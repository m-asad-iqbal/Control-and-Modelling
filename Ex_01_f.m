clc
clear
close all
format long
p = 800;
r1 = 2.6684;
r2 = 3.3620;
fi = 0.1592;
a = 996.2;


V1 = 4*pi*r1^3/3;
m1 = 800*V1;
V2 = 4*pi*(r2^3-r1^3)/3;
m2 = p*V2;
cp = 0.7;
A1 = 4*pi*r1^2;
gamma = 1e-10;
a11 = -a*A1/(m1*cp);
a12 = a*A1/(m1*cp);
b1 = gamma/(m1*cp);
a21 = a*A1/(m2*cp);
a22 = -a*A1/(m2*cp)-fi/m2;

A = [a11 a12; a21 a22];
B = [b1; 0];
C = [0 1];

[Num, Den] = tfdata(tf(ss(A,B,C,0)),'v');
c = Den(end);
Den = Den/c;
Num = Num/c;
tau = sqrt(Den(1));
zeta = Den(2)/(2*tau);
kp = Num(end);

G2 = tf(Num,Den);

K = kp;
tau_c = 7.8e5;
tau_d = 5e4;

KP = tau_c/(K*tau_d);
tau_i = tau_c;

s = tf('s');
PI = KP + 1/(tau_i*s);

Gcl_2 = feedback(G2*PI,1);
G1 = G2/gamma;

Gcl_1 = G1/(1+PI*G2);

[Num1, Den1] = tfdata(Gcl_1,'v');
[Num2, Den2] = tfdata(Gcl_2,'v');
syms s
G1_sym = poly2sym(Num1,s)/poly2sym(Den1,s);

T2_1s = -1e7*G1_sym/s; 
T2_1t = ilaplace(T2_1s)


G2_sym = poly2sym(Num2,s)/poly2sym(Den2,s);
T2_2s = G2_sym/s; 
T2_2t = ilaplace(T2_2s)

