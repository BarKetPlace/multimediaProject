clear all
close all
clc

N = 1000;
f = 10;
t = 2*pi*[0:1/N:1-1/N];
x  = sin(2*pi*f*t);
varn = .5;
n = sqrt(varn)*randn(1,N);

y = x + n;

%Probability that an sample belongs to the clean signal
q = .5;
    

figure, plot(y);
