#!/bin/bash -ex

# title notice
echo "HTK training for TIMIT from Cantab Research"
echo "See http://www.cantabResearch.com/HTKtimit.html"

# (C) Cantab Research: Permission is granted for free distribution and
# usage provided the copyright notice and the title notice are not
# altered.
#
# version 1.0: 22Jun06 - first public release
# version 1.1: 03Jul06 - vFloors was calculated but not used (Andrew Morris)
#                      - changed name from MODEL to MMF
#                      - fixed spurious HERest WARNING [-2331]
# version 1.2: 23Aug06 - fixed insertion errors caused by h# in references
# version 1.3: 14Sep06 - fixed typo KFLMAP -> KFLCFG in HResults line
#
# based on the Matlab scipts by Dan Ellis

# you *must* edit this block to reflect your local setup
TIMIT=/home/sach/HTK/TIMIT
SAMPLES=/home/sach/HTK/HTK_TIMIT/HTK_TIMIT1/samples

# do all work in this sub directory
WORK=work`echo $* | tr -d ' '`
mkdir -p $WORK

# uncomment this line to log everything in the working directory
# exec >& $WORK/$0.log

# some options to play with once the base system is working
HTKMFCC=true             # use HTK's MFCC instead of external USER executable
NMIXMONO=20              # maximum number of Gaussians per state in monophones
NMIXTRI=20               # maximum number of Gaussians per state in triphones
MINTESTMONO=1            # test the monophones after this number of Gaussians
MINTESTTRI=1             # test the triophones after this number of Gaussians
NPASSPERMIX=4            # number of fwd/bwd passes per mixture increase
TESTSET=coreTest         # set to "test" for full test set or "coreTest"
KFLMAP=false             # set to true to addionally output KFL mapped scores
echo "Start at `date`"

## keep a copy of this script,  plp source and executable
cp -p $0  $WORK
if ! $HTKMFCC ; then
  cp -p /cantab/src/plp.c `which plp` $WORK
fi

cd $WORK

# write the timit config used if using HTK MFCC
if $HTKMFCC ; then
  cat <<"EOF" > config
SOURCEKIND     = WAVEFORM
SOURCEFORMAT   = NIST
SAVECOMPRESSED = TRUE
SAVEWITHCRC    = TRUE
TARGETKIND     = MFCC_E_D_A_Z
TARGETRATE     = 100000
SOURCERATE     = 625
WINDOWSIZE     = 250000.0
PREEMCOEF      = 0.97
ZMEANSOURCE    = TRUE
USEHAMMING     = TRUE
CEPLIFTER      = 22
NUMCHANS       = 26
NUMCEPS        = 12
ENORMALISE     = TRUE
ESCALE         = 1.0
EOF
else
  cat <<"EOF" > config
HPARM: SOURCEFORMAT = HTK   
# HPARM: SOURCEKIND = USER
# HPARM: TARGETKIND = USER_D_A
EOF
fi

# read the TIMIT disk and encode into acoutic features
for DIR in train test ; do
  # create a mirror of the TIMIT directory structure
  (cd $TIMIT ; find -iname $DIR -type d) | xargs mkdir -p

  # generate lists of files
  (cd $TIMIT ; find -iname $DIR -type f -name s[ix]\*wav) | sort > $DIR.wav
  sed "s/wav$/phn/" $DIR.wav > $DIR.phn
  sed "s/wav$/mfc/" $DIR.wav > $DIR.scp
  sed "s/wav$/txt/" $DIR.wav > $DIR.txt

  # generate the acoutic feature vectors
  if $HTKMFCC ; then
    paste $DIR.wav $DIR.scp | sed "s:^:$TIMIT/:" > $DIR.convert
    HCopy -C config -S $DIR.convert
    rm -f $DIR.convert
  else
    sed "s/.wav$//" $DIR.wav | while read base ; do
      wavAddHead -skip 1024 $TIMIT/$base.wav tmp.wav
      plp -outputHTK tmp.wav tmp.plp
      deltas $* tmp.plp $base.mfc
    done
    rm -f tmp.wav
  fi

  # ax-h conflicts with HTK's triphone naming convention, so change it
  # also generate .txt files suitable for use in language modelling
  sed "s/.wav$//" $DIR.wav | while read base ; do
    sed 's/ ax-h$/ axh/' < $TIMIT/$base.phn > $base.phn
    egrep -v 'h#$' $base.phn > $base.txt
  done

  # create MLF
  HLEd -S $DIR.phn -i ${DIR}Mono.mlf /dev/null

  rm -f $DIR.wav
done

# filter the main test set to get the core test set
FILTER='^test/dr./[mf](DAB0|WBT0|ELC0|TAS1|WEW0|PAS0|JMP0|LNT0|PKT0|LLL0|TLS0|JLM0|BPM0|KLT0|NLP0|CMJ0|JDH0|MGD0|GRT0|NJM0|DHC0|JLN0|PAM0|MLD0)/s[ix]'
egrep -i $FILTER test.scp > coreTest.scp
egrep -i $FILTER test.phn > coreTest.phn
HLEd -S coreTest.phn -i coreTestMono.mlf /dev/null

# create list of monophones
find -iname train -name \*phn | xargs cat | awk '{print $3}' | sort -u > monophones

# and derive the dictionary and the modified monophone list from the monophones
egrep -v h# monophones > monophones-h#
paste monophones-h# monophones-h# > dict
echo "!ENTER	[] h#" >> dict
echo "!EXIT	[] h#" >> dict
echo '!ENTER' >> monophones-h#
echo '!EXIT' >> monophones-h#


# generate a template for a prototype model
cat <<"EOF" > sim.pcf
<BEGINproto_config_file>
<COMMENT>
   This PCF produces a single mixture, single stream prototype system
<BEGINsys_setup>
hsKind: P
covKind: D
nStates: 3
nStreams: 1
sWidths: 39
mixes: 1
parmKind: MFCC_E_D_A_Z
vecSize: 39
outDir: .
hmmList: protolist
<ENDsys_setup>
<ENDproto_config_file>
EOF

if ! $HTKMFCC ; then
  sed 's/^parmKind: MFCC_E_D_A_Z$/parmKind: USER_D_A/' < sim.pcf > tmp.pcf
  mv tmp.pcf sim.pcf
fi

# generate a prototype model
echo proto > protolist
echo N | $SAMPLES/HTKDemo/MakeProtoHMMSet sim.pcf

if ! $HTKMFCC ; then
  CONFIG="-C config"
fi

HCompV $CONFIG  -m -S train.scp -f 0.01 -M . -o new proto

KFLCFG='-e n en -e aa ao -e ah ax-h -e ah ax -e ih ix -e l el -e sh zh -e uw ux -e er axr -e m em -e n nx -e ng eng -e hh hv -e pau pcl -e pau tcl -e pau kcl -e pau q -e pau bcl -e pau dcl -e pau gcl -e pau epi -e pau h#'

nmix=1
NEWDIR=mono-nmix$nmix-npass0

# concatenate prototype models to build a flat-start model
mkdir -p $NEWDIR
sed '1,3!d' new > $NEWDIR/MMF
cat vFloors >> $NEWDIR/MMF
for i in `cat monophones` ; do
  sed -e "1,3d" -e "s/new/$i/" new >> $NEWDIR/MMF
done

HLEd -S train.txt -i trainTxt.mlf /dev/null
HLStats  -b bigfn -o -I trainTxt.mlf -S train.txt monophones-h#
HBuild  -n bigfn monophones-h# outLatFile

# HVITE needs to be told to perform cross word context expansion
cat <<"EOF" > hvite.config
FORCECXTEXP = TRUE
ALLOWXWRDEXP = TRUE
EOF

OPT="$CONFIG  -m 0 -t 250 150 1000 -S train.scp"

echo Start training monophones at: `date`


while [ $nmix -le $NMIXMONO ] ; do

  ## NB the inner loop of both cases is duplicated - change both! 
  if [ $nmix -eq 1 ] ; then
    npass=1;
    while [ $npass -le $NPASSPERMIX ] ; do
      OLDDIR=$NEWDIR
      NEWDIR=mono-nmix$nmix-npass$npass
      mkdir -p $NEWDIR
      HERest $OPT -I trainMono.mlf -H $OLDDIR/MMF -M $NEWDIR monophones > $NEWDIR/LOG
      npass=$(($npass+1))
    done
    echo 'MU 2 {*.state[2-4].mix}' > tmp.hed
    nmix=2
  else
    NEWDIR=mono-nmix$nmix-npass0 
    mkdir -p $NEWDIR
    HHEd -H $OLDDIR/MMF -M $NEWDIR tmp.hed monophones
    npass=1
    while [ $npass -le $NPASSPERMIX ] ; do
      OLDDIR=$NEWDIR
      NEWDIR=mono-nmix$nmix-npass$npass
      mkdir -p $NEWDIR
      HERest $OPT -I trainMono.mlf -H $OLDDIR/MMF -M $NEWDIR monophones > $NEWDIR/LOG
      npass=$(($npass+1))
    done
    echo 'MU +2 {*.state[2-4].mix}' > tmp.hed
    nmix=$(($nmix+2))
  fi

  # test models
  if [ $nmix -ge $MINTESTMONO ] ; then
    HVite $CONFIG -t 100 100 4000  -H $NEWDIR/MMF -S $TESTSET.scp -i $NEWDIR/recout.mlf -w outLatFile -p 0.0 -s 5.0 dict monophones
    HResults  -e '???' h# -I ${TESTSET}Mono.mlf monophones $NEWDIR/recout.mlf
    if $KFLMAP ; then
      HResults  -e '???' h# $KFLCFG -I ${TESTSET}Mono.mlf monophones $NEWDIR/recout.mlf >>RESULT_MONO
    fi
  fi

done

echo Completed monophone training at: `date`
