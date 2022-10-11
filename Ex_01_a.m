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

a11 = -a*A1/(m1*cp);
a12 = a*A1/(m1*cp);
b1 = 1/(m1*cp);
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

G1 = tf(Num,Den)
