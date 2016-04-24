#!/bin/bash

find /home/antoine/projectmultimedia/TIMIT/TRAIN -type f -iname *.WAV > train_list.txt
find /home/antoine/projectmultimedia/TIMIT/TEST -type f -iname *.WAV > test_list.txt

