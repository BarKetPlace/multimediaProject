function [] = funct_SpeechTimit(pathes,opt)
% Go build a structure of raw speech signals to be used in the remainder of the system
% Inputs :
%
% path train : path to a txt file that lists all training files
%
% path test : path to a txt file that lists all testing files
%
% opt.nbtrain : number of speech signals to convert / default :
%
%all
%
% opt.nbtest : number of speech signals to convert / default :
%
%all
%
% Outputs :
%
% Save the structures in a .mat in the current folder

if nargin<2
opt.nbtrain = 5000;
opt.nbtest = 5000;
end
% all speech signals to be converted
%Initialization
DATA = struct('utt',[],'rawSpeech',[],'frames',[],'mfcc',[],'part1',[],'ourFeature',[]);

fidTrain = fopen(pathes{1});
fidTest = fopen(pathes{2});
fidDev = fopen(pathes{3});
tlineTest = fgetl(fidTest);
tlineTrain = fgetl(fidTrain);
tlineDev = fgetl(fidDev);
%% TRAINING : Framing
i = 0;
while ischar(tlineTrain) && i < opt.nbtrain
i = i + 1;
%
% counter
% Name of utterance (for Kaldi tool)
C = strsplit(tlineTrain,{'/','.'}) ;
lenC = length(C);
% dataTrain.utt{i} = [C{9} ' ' C{10} ' ' C{11} ' ' C{12}] ;
% Name with train drX is blocking in Kaldi
% Trial without it:
DATA.utt{i} = [C{lenC-2} '_' C{lenC-1}] ;
% dataTrain.utt{i} = tlineTrain;
% Read Timit
y =readsph(tlineTrain,'s',-1);
% Add to structure
DATA.rawSpeech{i} = y' ;
% raw speech in row vector
% Update line
tlineTrain = fgetl(fidTrain);
end
save('dataTrain.mat','DATA','-v7.3') ;
clear DATA
%% TESTING : Framing
DATA = struct('utt',[],'rawSpeech',[],'frames',[],'mfcc',[],'part1',[],'ourFeature',[]);

i = 0;
% counter
while ischar(tlineTest) && i < opt.nbtest
i = i + 1;
% Name of utterance (for Kaldi tool)
C = strsplit(tlineTest,{'/','.'}) ;
lenC = length(C);
% dataTest.utt{i} = [C{9} ' ' C{10} ' ' C{11} ' ' C{12}] ;
% Name with train drX is blocking in Kaldi
% Trial without it:
% dataTest.utt{i} = [C{11} ' ' C{12}] ;
DATA.utt{i} = [C{lenC-2} '_' C{lenC-1}] ;
% dataTest.utt{i} = tlineTest;
% Read Timit
y =readsph(tlineTest,'s',-1);
% Add to structure
DATA.rawSpeech{i} = y' ;
% Raw speech in row vector
% Update line
tlineTest = fgetl(fidTest);
end
save('dataTest.mat','DATA','-v7.3') ;
clear DATA;
%% DEV : Framing
DATA = struct('utt',[],'rawSpeech',[],'frames',[],'mfcc',[],'part1',[],'ourFeature',[]);
i = 0;
while ischar(tlineDev)
i = i + 1;
%
% counter
% Name of utterance (for Kaldi tool)
C = strsplit(tlineDev,{'/','.'}) ;
lenC = length(C);
% dataDev.utt{i} = [C{9} ' ' C{10} ' ' C{11} ' ' C{12}] ;
% Name with train drX is blocking in Kaldi
% Trial without it:
% dataDev.utt{i} = [C{11} ' ' C{12}] ;
DATA.utt{i} = [C{lenC-2} '_' C{lenC-1}] ;
% dataDev.utt{i} = tlineDev;
% Read Timit
y =readsph(tlineDev,'s',-1);
% Add to structure
DATA.rawSpeech{i} = y' ;
% raw speech in row vector
% Update line
tlineDev = fgetl(fidDev);
end
save('dataDev.mat','DATA','-v7.3') ;
%% SAVING DATA
fclose(fidTrain) ;
fclose(fidTest) ;
fclose(fidDev) ;
end
