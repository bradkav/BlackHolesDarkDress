#!/bin/bash

BASEDIR=/home/kavanagh/PBH/
TARGETDIR=/home/kavanagh/PBH/sims/$1

RUNDIR=$TARGETDIR/run
OUTDIR=$TARGETDIR/out

echo "Setting up folder: $TARGETDIR"

mkdir -p $TARGETDIR
mkdir -p $RUNDIR
mkdir -p $OUTDIR

cp $BASEDIR/run/Gadget2 $RUNDIR

if [ "$2" = "lowres" ]
then
    cp $BASEDIR/run/PBH_lowres.param $RUNDIR/PBH.param
    echo "Using low-res settings..."
else
    cp $BASEDIR/run/PBH_hires.param $RUNDIR/PBH.param
    echo "Using high-res settings..."
fi

cp $BASEDIR/run/PBH1.dat $RUNDIR
cp $BASEDIR/run/ICs.txt $RUNDIR
cp $BASEDIR/RunPBH_template.sh $TARGETDIR/RunPBH_$1.sh

cd $TARGETDIR
#sed -i -- "s@RUN_DIRECTORY@$RUNDIR@g" RunPBH_$1.sh
sed -i -- "s@JOB_ID@$1@g" RunPBH_$1.sh

#cd $RUNDIR
#sed -i -- "s@RUN_DIRECTORY@$RUNDIR@g" PBH.param
#sed -i -- "s@OUT_DIRECTORY@$OUTDIR@g" PBH.param
 