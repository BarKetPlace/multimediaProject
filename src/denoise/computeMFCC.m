function computeMFCC(data_filename,snr,denoise_flag)
%function computeMFCC(data_filename,snr)
%Usage :: computeMFCC("dataTrain"/"dataTest"/"dataDev/dataTestNoisy5/10/15dB",snr)
%INPUT: data_filename   is the name of the datafile to process
%       snr             is the SNR of the noisy signal we want to process(put -1 if you want to process the clean data)
%       denoise_flag    indicates wether you want to use the denoising
%       algorithm or not
%OUTOUT: none


if nargin < 3
    denoise_flag=0;
end
%We put the ark files into the kaldi folder
    outputFolder = '/home/antoine/kaldi-trunk/egs/timit/s5/MatlabMFCC/';

if (strcmp(data_filename,'dataTrain'))%Compute MFCC on training datas
    ark_filename = [outputFolder 'raw_MatlabMFCC_train.ark'];

    
elseif (strcmp(data_filename,'dataTest'))%Compute MFCC on testing datas
    ark_filename = [outputFolder 'raw_MatlabMFCC_test.ark'];
        if (snr~=-1) %Compute on noisy features
            data_filename = [data_filename 'Noisy' num2str(snr) 'dB'];           
        end
        
elseif (strcmp(data_filename,'dataDev'))%Compute MFCC on dev datas
    ark_filename = [outputFolder 'raw_MatlabMFCC_dev.ark'];
    

else
    fprintf('Usage :: computeMFCC("dataTrain"/"dataTest"/"dataDev/dataTestNoisy5/10/15dB")\n');
    return;
end
    
    
%Then do the routine
% % Clean old MFCCs
if ~system(['test -f ' data_filename 'MFCC.mat']) %If file exists
    system(['rm ' data_filename 'MFCC.mat']);%Remove it
end
load([data_filename '.mat']);%Load data
if denoise_flag
    fprintf([data_filename ' MFCC Extraction & denoising...']);
else 
    fprintf([data_filename ' MFCC Extraction...']);
end

MFCCcell = getMFCC(DATA,denoise_flag);%extract MFCCs


fprintf('done.\n');
%dataTrain is very heavy so we create a copy containing only the utterance name
%and the MFCCs
DATA_MFCC = struct('utt',[],'rawSpeech',[],'frames',[],'feature',[],'part1',[],'ourFeature',[]);
DATA_MFCC.utt = DATA.utt;
clear DATA; %Useless
DATA_MFCC.feature = MFCCcell;
save([data_filename 'MFCC.mat'],'DATA_MFCC');
%   Convert & save to kaldi format
%   Clean the old ones if exists
if ~system(['test -f ' ark_filename]) %If file exists
    system(['rm ' ark_filename]);%Remove it
    system(['rm ' ark_filename(1:end-3) 'scp']);%Remove it
end

writekaldifeatures(DATA_MFCC,ark_filename);


end