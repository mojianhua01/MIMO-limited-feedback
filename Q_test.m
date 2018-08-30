clc;
clear;
close all;

%% 2-bit quantization

A = 6;
x_vector = 0.1:0.1:10
for ite = 1:1:length(x_vector)
    x = x_vector(ite);
    p1 = qfunc(A+x);
    p2 = qfunc(x) - qfunc(A+x);
    p3 = qfunc(x-A) - qfunc(x);
    p4 = 1 - qfunc(x-A);
   entropy(ite) = -p1 *log2(p1) - p2*log2(p2) - p3*log2(p3) - p4*log2(p4);
   [p1 p2 p3 p4]
end;

figure,
plot(x_vector, entropy)
grid on;
title('entropy')

figure,
plot(x_vector, -entropy)
grid on;
title('')