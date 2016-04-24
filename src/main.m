% clear all
close all
clc
clear all

%% Load database

% pathes{1} = '../train_list.txt';
% pathes{2} = '../test_list.txt';
% pathes{3} = '../test_list.txt';
% 
% funct_SpeechTimit(pathes);
% load dataTrain.mat
%% Load all files 
load dataTrain.mat
%% and process them
MFCCcell = getMFCC(dataTrain, FilterBank, frameLength);
dataTrainMFCC = struct('utt',[],'rawSpeech',[],'frames',[],'feature',[],'part1',[],'ourFeature',[]);
dataTrainMFCC.utt = dataTrain.utt;
clear dataTrain;
dataTrainMFCC.feature = MFCCcell;
save('dataTrainMFCC.mat','dataTrainMFCC','-v7.3');
%% Convert to kaldi format
load dataTrainMFCC.mat;
outputFolder = '/home/antoine/kaldi-trunk/egs/timit/s5/MatlabMFCC/';
filename = [outputFolder 'raw_MatlabMFCC_train.ark'];
%%
writekaldifeatures(dataTrainMFCC,filename)

% funct_WriteKaldiFormat(dataTrainMFCC.utt, dataTrainMFCC.mfcc, outputFolder,filename,NbFiles);
%% Plot MFCC
figure,
for i=1:nframe
    plot((i-1)*13+1:i*13,MFCC(i,:)); hold on;
end
title('frame by frame MFCCs');