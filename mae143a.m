close all
clear all
clc
%% Problem 1
%1a
disp('Setting up the Diophantine problem (solve a*x1+b*y1=f1 for x and y) in Problem 1a)')
b=[1 0 -20 0 64];
a=[1 0 -35 0 259 0 -225];
f1=[1 18 127 444 799 690 225];
[x1,y1] = RR_Diophantine(a,b,f1);
polesf1 = roots(x1);
zeroesf1 = roots(y1);

fprintf('The coefficients for x_1(s) is: [');
fprintf('%g, ', x1(1:end-1));
fprintf('%g]\n', x1(end));
fprintf('The coefficients for y_1(s) is: [ ');
fprintf('%g, ', y1(1:end-1));
fprintf('%g]\n', y1(end));
fprintf('The zeroes for D_1(s) is z = [');
fprintf('%g, ', zeroesf1(1:end-1));
fprintf('%g]\n', zeroesf1(end));
fprintf('The poles for D_1(s) is p = [');
fprintf('%g, ', polesf1(1:end-1));
fprintf('%g]\n', polesf1(end));
fprintf('\n')

%1b
disp('Setting up the Diophantine problem (solve a*x1+b*y1=f1 for x and y) in Problem 1b)')
b=[1 0 -20 0 64];
a=[1 0 -35 0 259 0 -225];
f2 = RR_PolyConv([1 1], [1 1], [1 3], [1 3], [1 5], [1 5]);
[x2,y2] = RR_Diophantine(a,b,f2);
p = 0;
while length(x2)<length(y2)
    f2 = conv(f2, [1 50]);
    p = p+1;
    [x2,y2]=RR_Diophantine(a,b,f2);
end
polesf2 = roots(x2);
zeroesf2= roots(y2);


fprintf('The minimum power, p, required to achieve a strictly proper system is: ')
fprintf('%g\n', p)
fprintf('The coefficients for x_2(s) is: [');
fprintf('%g, ', x2(1:end-1));
fprintf('%g]\n', x2(end));
fprintf('The coefficients for y_2(s) is: [ ');
fprintf('%g, ', y2(1:end-1));
fprintf('%g]\n', y2(end));
fprintf('The zeroes for D_2(s) is z = [');
fprintf('%g, ', zeroesf2(1:end-1));
fprintf('%g]\n', zeroesf2(end));
fprintf('The poles for D_2(s) is p = [');
fprintf('%g, ', polesf2(1:end-1));
fprintf('%g]\n', polesf2(end));
fprintf('\n')

%% Problem 3
s = tf('s');
d = 0.1;
sys = exp(-d*s);
n1 = 1;
n2 = 2;
n3 = 4;
n4 = 8;
%First Order Pade Approximation
f1s= pade(sys, n1)
%Second Order Pade Approximation
f2s= pade(sys, n2)
%Third Order Pade Approximation
f4s= pade(sys, n3)
%Fourth Order Pade Approximation
f8s= pade(sys, n4)

%eroes and Poles of System F1s
figure(1)
pzmap(f1s)
axis([-25 25 -1 1])
title('Pole-Zero Map for F_1(s)')
% Zeroes and Poles of System F2s
figure(2)
pzmap(f2s)
axis([-35 35 -20 20])
title('Pole-Zero Map for F_2(s)')
% Zeroes and Poles of System F4s
figure(3)
pzmap(f4s)
axis([-65 65 -60 60])
title('Pole-Zero Map for F_4(s)')
% Zeroes and Poles of System F8s
figure(4)
pzmap(f8s)
title('Pole-Zero Map for F_8(s)')

%% Problem 4
%Part a
figure(5)
rlocus(tf([1],[1 0]))
title('Root Locus for G(s)')

%Part 4b
%F_1(s)
figure(6)
rlocus(tf(f1s/s, [1 0]))
title('Root Locus For G_1(s), d = 0.1')
%F_2(s)
figure(7)
rlocus(tf(f2s/s, [1 0]))
title('Root Locus For G_2(s), d = 0.1')
%F_4(s)
figure(8)
rlocus(tf(f4s/s, [1 0]))
title('Root Locus For G_4(s), d = 0.1')
%F_8(s)
figure(9)
rlocus(tf(f8s/s, [1 0]))
title('Root Locus For G_8(s), d = 0.1')

%Part 4c 
%Using d =0.2 for all Pade Approximations
d = 0.2;
sys = exp(-d*s)
f1s= pade(sys, n1)

f2s= pade(sys, n2)

f4s= pade(sys, n3)

f8s= pade(sys, n4)

%F_1(s)
figure(10)
rlocus(tf(f1s/s, [1 0]))
title('Root Locus For G_1(s), d = 0.2')
%F_2(s)
figure(11)
rlocus(tf(f2s/s, [1 0]))
title('Root Locus For G_2(s), d = 0.2')
%F_4(s)
figure(12)
rlocus(tf(f4s/s, [1 0]))
title('Root Locus For G_4(s), d = 0.2')
%F_8(s)
figure(13)
rlocus(tf(f8s/s, [1 0]))
title('Root Locus For G_8(s), d = 0.2')
