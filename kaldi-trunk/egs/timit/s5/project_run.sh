#!/bin/bash

#
# Copyright 2013 Bagher BabaAli,
#           2014 Brno University of Technology (Author: Karel Vesely)
#
# TIMIT, description of the database:
# http://perso.limsi.fr/lamel/TIMIT_NISTIR4930.pdf
#
# Hon and Lee paper on TIMIT, 1988, introduces mapping to 48 training phonemes, 
# then re-mapping to 39 phonemes for scoring:
# http://repository.cmu.edu/cgi/viewcontent.cgi?article=2768&context=compsci



#To run this script, go into the kaldi/egs/timit/s5 directory
#
#Usage: Simplest: ./project_run.sh
#Using options:   ./project_run.sh --training true(default)/false --snr 5 or 10 or 15
#
#The variable training is used to skip :
#	the mfcc features extraction of the training data;
#	The training of the HMMs
#
#The variable snr is used to specify the SNR in dB of the test data. If you want to run kaldi on
#the clean data, just omit the option --snr

start=`date +%s` #Starting Time
startday=`date +%c` #Starting date
#Customize these 2 variables
timit=/home/antoine/Documents/multimediaProject/TIMIT
matlabcode=/home/antoine/Documents/multimediaProject/src/
##########################################################
kaldi_s5=`pwd`


comment="Clean testing data" #if no option specified, use the clean data
training=true #By default, extract the MFCC from training data and train the hmm
snr=-1 #Default, no noise on the testing data
denoise=false

. parse_options.sh || exit 1;

#Create the list of train, test and dev files. Store into data_[...].list in $matlabcode
if ! [ -f $matlabcode/dataTrain.mat ];then

cd $matlabcode; ./listfiles.sh $timit $kaldi_s5; cd $kaldi_s5;

#If the files dataTrain.mat, dataTest.mat and dataDev.mat do not exist in the matlab folder, run the following matlab script
echo "Load databases in files"
matlab -nojvm -nodesktop -r "cd('$matlabcode');\
				loadDB('data_train.list');\
				loadDB('data_test.list');\
				loadDB('data_dev.list');\
				quit;\
				"
fi

if [ $snr != -1 ];then
	comment="Noisy data, snr= $snr"
	#create the noisy data
	matlab -nojvm -nodesktop -r "cd('$matlabcode'); makenoise($snr); quit;"
fi

if [ $denoise == true ]; then
	comment=$comment" with denoising"
else
	comment=$comment" without denoising"
fi
. ./cmd.sh 
[ -f path.sh ] && . ./path.sh
set -e

# Acoustic model parameters
numLeavesTri1=2500
numGaussTri1=15000
numLeavesMLLT=2500
numGaussMLLT=15000
numLeavesSAT=2500
numGaussSAT=15000
numGaussUBM=400
numLeavesSGMM=7000
numGaussSGMM=9000

feats_nj=10
train_nj=30
decode_nj=10

echo ============================================================================
echo "                Data & Lexicon & Language Preparation                     "
echo ============================================================================
local/timit_data_prep.sh $timit || exit 1
local/timit_prepare_dict.sh

# Caution below: we remove optional silence by setting "--sil-prob 0.0",
# in TIMIT the silence appears also as a word in the dictionary and is scored.
utils/prepare_lang.sh --sil-prob 0.0 --position-dependent-phones false --num-sil-states 3 \
 data/local/dict "sil" data/local/lang_tmp data/lang

local/timit_format_data.sh

			    

mfccdir=MatlabMFCC

if [ ${training} == true ]; then
echo ============================================================================
echo "         Training data MFCC Feature Extration                             "
echo ============================================================================

#Call matlab to compute the mfcc on the training data
x=train
matlab -nojvm -nodesktop -r "cd('$matlabcode'); computeMFCC('dataTrain'); quit;"

#prepare the folders to match the kaldi requirements
sort "$mfccdir"/raw_${mfccdir}_$x.scp > data/feats.scp
cp data/feats.scp data/$x/feats.scp
mkdir -p exp/make_$mfccdir
mkdir -p exp/mono

steps/compute_cmvn_stats.sh data/$x exp/make_$mfccdir/$x ./$mfccdir
utils/validate_data_dir.sh data/$x/


echo ============================================================================
echo "         Monophone training       "
echo ============================================================================
#Perform the monophone training
steps/train_mono.sh  --nj "$train_nj" --cmd "$train_cmd" data/train data/lang exp/mono
utils/mkgraph.sh --mono data/lang_test_bg exp/mono exp/mono/graph
fi


echo ============================================================================
echo "         Testing data MFCC Feature Extration                             "
echo ============================================================================
x=test

#   steps/make_mfcc.sh --cmd "$train_cmd" --nj $feats_nj data/$x exp/make_mfcc/$x $mfccdir
#Compute MFCC for the testing batch
matlab -nojvm -nodesktop -r "cd('$matlabcode'); computeMFCC('dataTest',$snr,$denoise); quit;"

#prepare folders for kaldi requirements
sort "$mfccdir"/raw_${mfccdir}_$x.scp > data/feats.scp
cp data/feats.scp data/$x/feats.scp
mkdir -p exp/make_$mfccdir
mkdir -p exp/mono

steps/compute_cmvn_stats.sh data/$x exp/make_$mfccdir/$x ./$mfccdir
utils/validate_data_dir.sh data/$x/


echo ============================================================================
echo "                     mono Phone Decoding                        "
echo ============================================================================



#steps/decode.sh --nj "$decode_nj" --cmd "$decode_cmd" \
# exp/mono/graph data/dev exp/mono/decode_dev

steps/decode.sh --nj "$decode_nj" --cmd "$decode_cmd" \
 exp/mono/graph data/test exp/mono/decode_test

#echo ============================================================================
#echo "           tri1 : Deltas + Delta-Deltas Training & Decoding               "
#echo ============================================================================

#steps/align_si.sh --boost-silence 1.25 --nj "$train_nj" --cmd "$train_cmd" \
# data/train data/lang exp/mono exp/mono_ali

# Train tri1, which is deltas + delta-deltas, on train data.
#steps/train_deltas.sh --cmd "$train_cmd" \
# $numLeavesTri1 $numGaussTri1 data/train data/lang exp/mono_ali exp/tri1

#utils/mkgraph.sh data/lang_test_bg exp/tri1 exp/tri1/graph

#steps/decode.sh --nj "$decode_nj" --cmd "$decode_cmd" \
# exp/tri1/graph data/dev exp/tri1/decode_dev

#steps/decode.sh --nj "$decode_nj" --cmd "$decode_cmd" \
# exp/tri1/graph data/test exp/tri1/decode_test

#echo ============================================================================
#echo "                 tri2 : LDA + MLLT Training & Decoding                    "
#echo ============================================================================

#steps/align_si.sh --nj "$train_nj" --cmd "$train_cmd" \
#  data/train data/lang exp/tri1 exp/tri1_ali

#steps/train_lda_mllt.sh --cmd "$train_cmd" \
# --splice-opts "--left-context=3 --right-context=3" \
# $numLeavesMLLT $numGaussMLLT data/train data/lang exp/tri1_ali exp/tri2

#utils/mkgraph.sh data/lang_test_bg exp/tri2 exp/tri2/graph

#steps/decode.sh --nj "$decode_nj" --cmd "$decode_cmd" \
# exp/tri2/graph data/dev exp/tri2/decode_dev

#steps/decode.sh --nj "$decode_nj" --cmd "$decode_cmd" \
# exp/tri2/graph data/test exp/tri2/decode_test

#echo ============================================================================
#echo "              tri3 : LDA + MLLT + SAT Training & Decoding                 "
#echo ============================================================================

## Align tri2 system with train data.
#steps/align_si.sh --nj "$train_nj" --cmd "$train_cmd" \
# --use-graphs true data/train data/lang exp/tri2 exp/tri2_ali

## From tri2 system, train tri3 which is LDA + MLLT + SAT.
#steps/train_sat.sh --cmd "$train_cmd" \
# $numLeavesSAT $numGaussSAT data/train data/lang exp/tri2_ali exp/tri3

#utils/mkgraph.sh data/lang_test_bg exp/tri3 exp/tri3/graph

#steps/decode_fmllr.sh --nj "$decode_nj" --cmd "$decode_cmd" \
# exp/tri3/graph data/dev exp/tri3/decode_dev

#steps/decode_fmllr.sh --nj "$decode_nj" --cmd "$decode_cmd" \
# exp/tri3/graph data/test exp/tri3/decode_test

#echo ============================================================================
#echo "                        SGMM2 Training & Decoding                         "
#echo ============================================================================

#steps/align_fmllr.sh --nj "$train_nj" --cmd "$train_cmd" \
# data/train data/lang exp/tri3 exp/tri3_ali

##exit 0 # From this point you can run Karel's DNN : local/nnet/run_dnn.sh 

#steps/train_ubm.sh --cmd "$train_cmd" \
# $numGaussUBM data/train data/lang exp/tri3_ali exp/ubm4

#steps/train_sgmm2.sh --cmd "$train_cmd" $numLeavesSGMM $numGaussSGMM \
# data/train data/lang exp/tri3_ali exp/ubm4/final.ubm exp/sgmm2_4

#utils/mkgraph.sh data/lang_test_bg exp/sgmm2_4 exp/sgmm2_4/graph

#steps/decode_sgmm2.sh --nj "$decode_nj" --cmd "$decode_cmd"\
# --transform-dir exp/tri3/decode_dev exp/sgmm2_4/graph data/dev \
# exp/sgmm2_4/decode_dev

#steps/decode_sgmm2.sh --nj "$decode_nj" --cmd "$decode_cmd"\
# --transform-dir exp/tri3/decode_test exp/sgmm2_4/graph data/test \
# exp/sgmm2_4/decode_test

#echo ============================================================================
#echo "                    MMI + SGMM2 Training & Decoding                       "
#echo ============================================================================

#steps/align_sgmm2.sh --nj "$train_nj" --cmd "$train_cmd" \
# --transform-dir exp/tri3_ali --use-graphs true --use-gselect true \
# data/train data/lang exp/sgmm2_4 exp/sgmm2_4_ali

#steps/make_denlats_sgmm2.sh --nj "$train_nj" --sub-split "$train_nj" \
# --acwt 0.2 --lattice-beam 10.0 --beam 18.0 \
# --cmd "$decode_cmd" --transform-dir exp/tri3_ali \
# data/train data/lang exp/sgmm2_4_ali exp/sgmm2_4_denlats

#steps/train_mmi_sgmm2.sh --acwt 0.2 --cmd "$decode_cmd" \
# --transform-dir exp/tri3_ali --boost 0.1 --drop-frames true \
# data/train data/lang exp/sgmm2_4_ali exp/sgmm2_4_denlats exp/sgmm2_4_mmi_b0.1

#for iter in 1 2 3 4; do
#  steps/decode_sgmm2_rescore.sh --cmd "$decode_cmd" --iter $iter \
#   --transform-dir exp/tri3/decode_dev data/lang_test_bg data/dev \
#   exp/sgmm2_4/decode_dev exp/sgmm2_4_mmi_b0.1/decode_dev_it$iter

#  steps/decode_sgmm2_rescore.sh --cmd "$decode_cmd" --iter $iter \
#   --transform-dir exp/tri3/decode_test data/lang_test_bg data/test \
#   exp/sgmm2_4/decode_test exp/sgmm2_4_mmi_b0.1/decode_test_it$iter
#done

#echo ============================================================================
#echo "                    DNN Hybrid Training & Decoding                        "
#echo ============================================================================

## DNN hybrid system training parameters
#dnn_mem_reqs="--mem 1G"
#dnn_extra_opts="--num_epochs 20 --num-epochs-extra 10 --add-layers-period 1 --shrink-interval 3"

#steps/nnet2/train_tanh.sh --mix-up 5000 --initial-learning-rate 0.015 \
#  --final-learning-rate 0.002 --num-hidden-layers 2  \
#  --num-jobs-nnet "$train_nj" --cmd "$train_cmd" "${dnn_train_extra_opts[@]}" \
#  data/train data/lang exp/tri3_ali exp/tri4_nnet

#[ ! -d exp/tri4_nnet/decode_dev ] && mkdir -p exp/tri4_nnet/decode_dev
#decode_extra_opts=(--num-threads 6)
#steps/nnet2/decode.sh --cmd "$decode_cmd" --nj "$decode_nj" "${decode_extra_opts[@]}" \
#  --transform-dir exp/tri3/decode_dev exp/tri3/graph data/dev \
#  exp/tri4_nnet/decode_dev | tee exp/tri4_nnet/decode_dev/decode.log

#[ ! -d exp/tri4_nnet/decode_test ] && mkdir -p exp/tri4_nnet/decode_test
#steps/nnet2/decode.sh --cmd "$decode_cmd" --nj "$decode_nj" "${decode_extra_opts[@]}" \
#  --transform-dir exp/tri3/decode_test exp/tri3/graph data/test \
#  exp/tri4_nnet/decode_test | tee exp/tri4_nnet/decode_test/decode.log

#echo ============================================================================
#echo "                    System Combination (DNN+SGMM)                         "
#echo ============================================================================

#for iter in 1 2 3 4; do
#  local/score_combine.sh --cmd "$decode_cmd" \
#   data/dev data/lang_test_bg exp/tri4_nnet/decode_dev \
#   exp/sgmm2_4_mmi_b0.1/decode_dev_it$iter exp/combine_2/decode_dev_it$iter

#  local/score_combine.sh --cmd "$decode_cmd" \
#   data/test data/lang_test_bg exp/tri4_nnet/decode_test \
#   exp/sgmm2_4_mmi_b0.1/decode_test_it$iter exp/combine_2/decode_test_it$iter
#done

#echo ============================================================================
#echo "               DNN Hybrid Training & Decoding (Karel's recipe)            "
#echo ============================================================================

#local/nnet/run_dnn.sh
##local/nnet/run_autoencoder.sh : an example, not used to build any system,

echo ============================================================================
echo "                    Getting Results [see RESULTS file]                    "
echo ============================================================================
stop=`date +%s` #End time
secs=`expr $stop - $start` #Elapsed time

echo >> RESULTS #Skip a line in file 
echo "============================ Comments: " $comment>> RESULTS 
echo "Started on $startday "========================================================== >> RESULTS
bash RESULTS test >> RESULTS #Write in file
bash RESULTS test # write in command line
echo "===========================Finished in `printf '%dh:%dm:%ds\n' $(($secs/3600)) $(($secs%3600/60)) $(($secs%60))`=================================" >> RESULTS

exit 0
