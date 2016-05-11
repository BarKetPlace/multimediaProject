% This function is used to evaluate the features which is subsequently used
% by the HTK system

function noise_file(input,SNR,Noise,train_flag)







s=input;


s=s-mean(s);

s= s / max(abs(s));


% ------------ Corrupting by noise ----------------


Noise_For_Input=Noise(1:length(s));

if (train_flag == 1)  % Then traning and no input noise (clean condition training)
    
    s_noisy=s;
    
else
    
    UV_Noise_For_Input=Noise_For_Input/std(Noise_For_Input);
    
    noise_var=(var(s) / (10^(SNR/10)));
    Noise_sample= (noise_var^(0.5)) * UV_Noise_For_Input;
    
    s_noisy=s+Noise_sample;
    
end


wavwrite(s,16000,'clean.wav');

wavwrite(s_noisy,16000,'clean_babble_snr0.wav');

plot(s);
figure; plot(s_noisy);



