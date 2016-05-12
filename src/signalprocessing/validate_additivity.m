clear all
cd /home/antoine/Documents/multimediaProject/src/signalprocessing

%% 
isignal= 56;
SNR=5;
noise_path = '/home/antoine/Documents/multimediaProject/TIMIT/NoiseDB/NoiseX_16kHz/';
noise_file = 'white_16kHz.wav';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load dataTest.mat

x = DATA.rawSpeech{1,isignal};
clear DATA

[Noise, Fs] = audioread([noise_path noise_file]);

%Conditioning of x
x=x-mean(x);
x=x/max(abs(x));
%extract the right noise length

noise_sig = Noise(1:length(x))';
UVnoise_sig = noise_sig/std(noise_sig);
varn= (var(x)) / (10^(SNR/10));
noise_sig = (varn^(.5))*UVnoise_sig;
noise_sig = noise_sig-mean(noise_sig);
n=noise_sig;
y=x+n;


framelen=.025;
cd ..
[cepstrax,aspectrumx,pspectrumx] = melfcc(x, Fs, [],'useenergy',1);
[cepstray,aspectrumy,pspectrumy] = melfcc(y, Fs, [],'useenergy',1);
[cepstran,aspectrumn,pspectrumn] = melfcc(n, Fs, [],'useenergy',1);
cd signalprocessing

nbframe=size(cepstrax,2);


for iframe = 1:nbframe
Ey= aspectrumy(:,iframe);
Ex= aspectrumx(:,iframe);
En= aspectrumn(:,iframe);
snr_mel_energy(iframe)= 10*log10(var(Ex)/var(En));
timewindow=round(1+(iframe-1)*framelen*Fs:min(length(x),iframe*framelen*Fs));
snr_time(iframe)=10*log10(var(x(timewindow))/var(n(timewindow)));
end

 
%%

iframe=round( (nbframe-1)*(rand()) ) +1;
iframe=33;
Ey= aspectrumy(:,iframe);
Ex= aspectrumx(:,iframe);
En= aspectrumn(:,iframe);
logPy=exp(cepstray(1,:));
% logPx=logPx+abs(min(logPx));

threshold = min(logPy)+.28*abs(max(logPy) - min(logPy));

if logPy(iframe)<=threshold
    decision='Silence frame';
else
    decision='Speech frame';
end


figure(1), clf
plot(Ey,'LineWidth',2); hold on;
plot(Ex+En,'LineWidth',2); hold on
plot(Ex); hold on;
plot(En);
legend('Ey','Ex+En','Ex','En');
xlabel('Mel space coefficent');
ylabel('Value f coefficient');
title({['mel SNR, 10log10(var(Ex)/var(En)):: ' num2str(snr_mel_energy(iframe))];...
        ['time SNR, 10log10(var(xframe)/var(nframe)):: ' num2str(snr_time(iframe))]});

figure(2), clf
plot(snr_mel_energy); hold on;
plot(snr_time);

legend('mel energy domain','time domain');
title('Frame by frame SNR');
xlabel('Frame number');
ylabel('SNRdB');


figure(3), clf
plot(logPy); hold on;
plot(iframe*[1 1], [min(logPy) max(logPy)]);hold on;
% plot([1 nbframe], threshold*[1 1]);
xlabel('Frame number');
ylabel('sum(abs(X(f))^2)');
title('Power of each frame');

energy_snr=[snr_mel_energy;logPy];
%%
[~,idx] = sort(snr_mel_energy);


figure(4), clf;
plot(snr_mel_energy(idx),logPy(idx))
ylabel('Py power of a frame');
xlabel('SNR in mel domain');