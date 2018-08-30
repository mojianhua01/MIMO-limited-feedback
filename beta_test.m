clc;
clear all;
syms x
eta = 1;
sigma = 1;
m = 10;
int(1/(eta *x + sigma) *(m-1) * (1-x)^(m-2),x,0,1)