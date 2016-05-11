% This file is to create a noise x database of Fs=16 kHz and ZMUV.


clear;
close all;

Fs=8000;

FilterLength=100;

input=wavread('~/RESEARCH3/NoiseDB/NoiseX/volvo.wav');

s1=resample(input,32,37,FilterLength);

s2=resample(s1,25,27,FilterLength);

s2=resample(s1,1,2,FilterLength);

s2=s2(100000:300000);


if (max(abs(s2)) < 1)
    s3=s2;
else
    s3=s2 / max(abs(s2));
end

wavwrite(s3,8000,'volvo_8kHz.wav');




