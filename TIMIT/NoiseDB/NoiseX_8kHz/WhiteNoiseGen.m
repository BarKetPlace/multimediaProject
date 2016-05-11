% This file provides a noise sample which is ZMUV Gaussian

clear;
close all;

Fs=8000;

SampleLength=500000;

noise=randn(1,SampleLength);

noise1=noise(100000:300000);

if (max(abs(noise1)) < 1)
    noise2=noise1;
else
    noise2=noise1 / max(abs(noise1));
end

wavwrite(noise2,Fs,'white_8kHz.wav');
