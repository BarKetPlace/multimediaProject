#!/bin/bash
timit=$1; #/home/antoine/Documents/multimediaProject/TIMIT;
kaldi_egs=$2; #/home/antoine/kaldi-trunk/egs/timit/s5;

#training data list
find $timit/TRAIN -type f -iname *.WAV > train_list.txt

#	Do not take into account all the files (to match the kaldi/egs/timit/s5/run.sh script)
sed -n '/SA1/!p' train_list.txt > train_list_tmp.txt
sed -n '/SA2/!p' train_list_tmp.txt > data_train.list

#Testing datalist
find $timit/TEST -type f -iname *.WAV > test_list.txt

#	Do not take into account all the files (to match the kaldi/egs/timit/s5/run.sh script)
sed -n '/SA1/!p' test_list.txt > test_list_tmp.txt
sed -n '/SA2/!p' test_list_tmp.txt > data_test.list.tmp
#	Match the list of speakers from kaldi example
#	The list is conf/test_spk.list, we need to make it upper case
sed 's/\([a-z]\)/\U\1/g' $kaldi_egs/conf/test_spk.list > $kaldi_egs/conf/test_SPK.list
#	We match the speaker list with our test_list.txt
grep -f $kaldi_egs/conf/test_SPK.list data_test.list.tmp > data_test.list


#Devlopment datalist
cp test_list.txt dev_list.txt
#Same processing as previously
sed -n '/SA1/!p' dev_list.txt > dev_list_tmp.txt
sed -n '/SA2/!p' dev_list_tmp.txt > data_dev.list.tmp
#	Match the list of speakers from kaldi example
#	The list is conf/test_spk.list, we need to make it upper case
sed 's/\([a-z]\)/\U\1/g' $kaldi_egs/conf/dev_spk.list > $kaldi_egs/conf/dev_SPK.list
#	We match the speaker list with our dev_list.txt
grep -f $kaldi_egs/conf/dev_SPK.list data_dev.list.tmp > data_dev.list

#clean up
rm train_list*  test_list* dev_list* *.tmp
