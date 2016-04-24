% clear all
close all
clc
clear all

%% Load databases 
if (system('test -f dataTrain.mat')) %If file doesnot exists
    system('./listfiles.sh');%Create the list of files, should be modified according to system architecture
    pathes{1} = 'train_list.txt';%List of training files
    pathes{2} = 'test_list.txt';%List of testing files
    pathes{3} = 'test_list.txt';%List of dev files
    fprintf('Preparation of TIMIT database... ');
    funct_SpeechTimit(pathes);%Three .mat files are saved in the current folder
    fprintf('Done.\n');
end

%% Process training Data
fprintf('Loading Training data... ');
load dataTrain.mat;
fprintf('Done.\n');
fprintf('MFCC Extraction... ');
MFCCcell = getMFCC(dataTrain);%extract MFCCs
fprintf('Done.\n');
%dataTrain is very heavy so we create a copy containing only the utterance name
%and the MFCCs
dataTrainMFCC = struct('utt',[],'rawSpeech',[],'frames',[],'feature',[],'part1',[],'ourFeature',[]);
dataTrainMFCC.utt = dataTrain.utt;
clear dataTrain; %Useless
dataTrainMFCC.feature = MFCCcell;
save('dataTrainMFCC.mat','dataTrainMFCC','-v7.3');
% Convert to kaldi format
outputFolder = '/home/antoine/kaldi-trunk/egs/timit/s5/MatlabMFCC/';
filename = [outputFolder 'raw_MatlabMFCC_train.ark'];
writekaldifeatures(dataTrainMFCC,filename);
%% Process testing data
fprintf('Loading Testing data ...\n');
load dataTest.mat;
fprintf('\bDone.\n');
fprintf('MFCC Extraction... ');
MFCCcell = getMFCC(dataTest);%extract MFCCs
fprintf('Done.\n');
%dataTest is very heavy so we create a copy containing only the utterance name
%and the MFCCs
dataTestMFCC = struct('utt',[],'rawSpeech',[],'frames',[],'feature',[],'part1',[],'ourFeature',[]);
dataTestMFCC.utt = dataTest.utt;
clear dataTest;
dataTestMFCC.feature = MFCCcell;
save('dataTestMFCC.mat','dataTestMFCC','-v7.3');
% Convert to kaldi format
outputFolder = '/home/antoine/kaldi-trunk/egs/timit/s5/MatlabMFCC/';
filename = [outputFolder 'raw_MatlabMFCC_test.ark'];
writekaldifeatures(dataTestMFCC,filename);
%% Process dev data
fprintf('Loading Dev data... ');
load dataDev.mat;
fprintf('Done.\n');
fprintf('MFCC Extraction... ');
MFCCcell = getMFCC(dataDev);%extract MFCCs
fprintf('Done.\n');
%dataDev is very heavy so we create a copy containing only the utterance name
%and the MFCCs
dataDevMFCC = struct('utt',[],'rawSpeech',[],'frames',[],'feature',[],'part1',[],'ourFeature',[]);
dataDevMFCC.utt = dataDev.utt;
clear dataDev;
dataDevMFCC.feature = MFCCcell;
save('dataDevMFCC.mat','dataDevMFCC','-v7.3');
% Convert to kaldi format
outputFolder = '/home/antoine/kaldi-trunk/egs/timit/s5/MatlabMFCC/';
filename = [outputFolder 'raw_MatlabMFCC_Dev.ark'];
writekaldifeatures(dataDevMFCC,filename);