% This file is to create a noise x database of Fs=16 kHz and ZMUV.


clear;
close all;

Fs=16000;

FilterLength=100;

input=wavread('~/RESEARCH3/NoiseDB/NoiseX/hfchannel.wav');

s1=resample(input,32,37,FilterLength);

s2=resample(s1,25,27,FilterLength);

s2=s2(100000:300000);


if (max(abs(s2)) < 1)
    s3=s2;
else
    s3=s2 / max(abs(s2));
end

wavwrite(s3,16000,'hfchannel_16kHz.wav');




