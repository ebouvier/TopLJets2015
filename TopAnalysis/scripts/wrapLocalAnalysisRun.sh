#!/bin/bash

#determine CMSSW config
SCRIPT=$(readlink -f $0)
SCRIPTPATH=`dirname $SCRIPT`
ARCH=${SCRIPTPATH##/*/}
WORKDIR=${SCRIPTPATH}/../

#configure environment
cd $WORKDIR
export SCRAM_ARCH=$ARCH
eval `scram r -sh`

#run with the arguments passed
$*
#xrdcp MC13TeV_*.root root://eoscms//eos/cms/store/user/byates/LJets2015/skim/
#xrdcp Data13TeV_*.root root://eoscms//eos/cms/store/user/byates/LJets2015/skim/
