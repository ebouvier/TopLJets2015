#!/bin/bash

WHAT=$1; 
ERA=$2
if [ "$#" -ne 2 ]; then 
    echo "steerTOPMassAnalysis.sh <SEL/MERGESEL/PLOTSEL/WWWSEL> <ERA>";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGESEL     - merge the output of the jobs";
    echo "        PLOTSEL      - runs the plotter tool on the selection";
    echo "        WWWSEL       - moves the plots to an afs-web-based area";
    echo " "
    echo "        ERA          - era2015/era2016";
    exit 1; 
fi

export LSB_JOB_REPORT_MAIL=N

queue=8nh
githash=8db9ad6
lumi=12868.66
lumiUnc=0.062
#eosdir=/store/user/byates/LJets2015/${githash}
eosdir=/store/user/byates/LJets2015/8db9ad6/MC13TeV_TTJets_powheg/MergedMiniEvents_91.root

case $ERA in
    era2015)
	githash=8c1e7c9;
	lumi=2267.84
	eosdir=/store/cmst3/user/psilva/LJets2015/${githash}
	;;
esac

outdir=/afs/cern.ch/work/e/ebouvier/BMesons13TeV/CMSSW_8_0_11-Elvire/src/TopLJets2015/TopAnalysis/sel/
wwwdir=~/www/Elvire/


RED='\e[31m'
NC='\e[0m'
case $WHAT in
    SEL )
	python scripts/runLocalAnalysis.py -i ${eosdir} -n 8 -q ${queue} -o ${outdir} --era ${ERA} -m TOPMass::RunTopMass ;
	;;
    MERGESEL )
	./scripts/mergeOutputs.py ${outdir}/sel True;	
	;;
    PLOTSEL )
	python scripts/plotter.py -i ${outdir}/sel --puNormSF PU_WgtCtr_all  -j data/${ERA}/samples.json -l ${lumi} --saveLog;# --mcUnc ${lumiUnc};	
        #python scripts/plotter.py -i ${outdir} -j data/${ERA}/samples.json -l ${lumi};
	;;
    WWWSEL )
	cp -p  ${outdir}/sel/plots/*.{png,pdf} ${wwwdir}/
	mv ${wwwdir}/*_log.{png,pdf} ${wwwdir}/log/
	mv ${wwwdir}/*_ee*.{png,pdf} ${wwwdir}/ee/
	mv ${wwwdir}/*_em*.{png,pdf} ${wwwdir}/em/
	mv ${wwwdir}/*_e*.{png,pdf} ${wwwdir}/e/
	mv ${wwwdir}/*_mm*.{png,pdf} ${wwwdir}/mm/
	mv ${wwwdir}/*_m*.{png,pdf} ${wwwdir}/m/
	cp -p test/index.php ${wwwdir}/
	cp -p test/index.php ${wwwdir}/log/
	cp -p test/index.php ${wwwdir}/ee/
	cp -p test/index.php ${wwwdir}/e/
	cp -p test/index.php ${wwwdir}/em/
	cp -p test/index.php ${wwwdir}/mm/
	cp -p test/index.php ${wwwdir}/m/
	;;
esac
