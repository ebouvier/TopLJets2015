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
eosdir=/store/user/byates/LJets2015/${githash}
#eosdir=/store/user/byates/LJets2015/8db9ad6/MC13TeV_TTJets_powheg/MergedMiniEvents_91.root

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
	./scripts/mergeOutputs.py ${outdir} True;	
	;;
    PLOTSEL )
	python scripts/plotter.py -i ${outdir} --puNormSF puwgtctr  -j data/${ERA}/samples.json -l ${lumi} --saveLog;# --mcUnc ${lumiUnc};	
        #python scripts/plotter.py -i ${outdir} -j data/${ERA}/samples.json -l ${lumi};
	;;
    WWWSEL )
	cp -p  ${outdir}/plots/*.{png,pdf} ${wwwdir}/

	mv ${wwwdir}/*_check*.{png,pdf} ${wwwdir}/check/
	mv ${wwwdir}/check/*_no_weight*.{png,pdf} ${wwwdir}/check/unweighted/
	mv ${wwwdir}/check/unweighted/*_log*.{png,pdf} ${wwwdir}/check/unweighted/log/
	mv ${wwwdir}/check/unweighted/log/*_all*.{png,pdf} ${wwwdir}/check/unweighted/log/all/
	mv ${wwwdir}/check/unweighted/log/*_ee*.{png,pdf} ${wwwdir}/check/unweighted/log/ee/
	mv ${wwwdir}/check/unweighted/log/*_em*.{png,pdf} ${wwwdir}/check/unweighted/log/em/
	mv ${wwwdir}/check/unweighted/log/*_mm*.{png,pdf} ${wwwdir}/check/unweighted/log/mm/
	mv ${wwwdir}/check/unweighted/log/*_e*.{png,pdf} ${wwwdir}/check/unweighted/log/e/
	mv ${wwwdir}/check/unweighted/log/*_m*.{png,pdf} ${wwwdir}/check/unweighted/log/m/
	mv ${wwwdir}/check/unweighted/*.{png,pdf} ${wwwdir}/check/unweighted/norm/
	mv ${wwwdir}/check/unweighted/norm/*_all*.{png,pdf} ${wwwdir}/check/unweighted/norm/all/
	mv ${wwwdir}/check/unweighted/norm/*_ee*.{png,pdf} ${wwwdir}/check/unweighted/norm/ee/
	mv ${wwwdir}/check/unweighted/norm/*_em*.{png,pdf} ${wwwdir}/check/unweighted/norm/em/
	mv ${wwwdir}/check/unweighted/norm/*_mm*.{png,pdf} ${wwwdir}/check/unweighted/norm/mm/
	mv ${wwwdir}/check/unweighted/norm/*_e*.{png,pdf} ${wwwdir}/check/unweighted/norm/e/
	mv ${wwwdir}/check/unweighted/norm/*_m*.{png,pdf} ${wwwdir}/check/unweighted/norm/m/ 
  mv ${wwwdir}/check/*.{png,pdf} ${wwwdir}/check/weighted/
	mv ${wwwdir}/check/weighted/*_log*.{png,pdf} ${wwwdir}/check/weighted/log/
	mv ${wwwdir}/check/weighted/log/*_all*.{png,pdf} ${wwwdir}/check/weighted/log/all/
	mv ${wwwdir}/check/weighted/log/*_ee*.{png,pdf} ${wwwdir}/check/weighted/log/ee/
	mv ${wwwdir}/check/weighted/log/*_em*.{png,pdf} ${wwwdir}/check/weighted/log/em/
	mv ${wwwdir}/check/weighted/log/*_mm*.{png,pdf} ${wwwdir}/check/weighted/log/mm/
	mv ${wwwdir}/check/weighted/log/*_e*.{png,pdf} ${wwwdir}/check/weighted/log/e/
	mv ${wwwdir}/check/weighted/log/*_m*.{png,pdf} ${wwwdir}/check/weighted/log/m/
	mv ${wwwdir}/check/weighted/*.{png,pdf} ${wwwdir}/check/weighted/norm/
	mv ${wwwdir}/check/weighted/norm/*_all*.{png,pdf} ${wwwdir}/check/weighted/norm/all/
	mv ${wwwdir}/check/weighted/norm/*_ee*.{png,pdf} ${wwwdir}/check/weighted/norm/ee/
	mv ${wwwdir}/check/weighted/norm/*_em*.{png,pdf} ${wwwdir}/check/weighted/norm/em/
	mv ${wwwdir}/check/weighted/norm/*_mm*.{png,pdf} ${wwwdir}/check/weighted/norm/mm/
	mv ${wwwdir}/check/weighted/norm/*_e*.{png,pdf} ${wwwdir}/check/weighted/norm/e/
	mv ${wwwdir}/check/weighted/norm/*_m*.{png,pdf} ${wwwdir}/check/weighted/norm/m/ 
  
	mv ${wwwdir}/*_topsel*.{png,pdf} ${wwwdir}/topsel/
	mv ${wwwdir}/topsel/*_no_weight*.{png,pdf} ${wwwdir}/topsel/unweighted/
	mv ${wwwdir}/topsel/unweighted/*_log*.{png,pdf} ${wwwdir}/topsel/unweighted/log/
	mv ${wwwdir}/topsel/unweighted/log/*_all*.{png,pdf} ${wwwdir}/topsel/unweighted/log/all/
	mv ${wwwdir}/topsel/unweighted/log/*_ee*.{png,pdf} ${wwwdir}/topsel/unweighted/log/ee/
	mv ${wwwdir}/topsel/unweighted/log/*_em*.{png,pdf} ${wwwdir}/topsel/unweighted/log/em/
	mv ${wwwdir}/topsel/unweighted/log/*_mm*.{png,pdf} ${wwwdir}/topsel/unweighted/log/mm/
	mv ${wwwdir}/topsel/unweighted/log/*_e*.{png,pdf} ${wwwdir}/topsel/unweighted/log/e/
	mv ${wwwdir}/topsel/unweighted/log/*_m*.{png,pdf} ${wwwdir}/topsel/unweighted/log/m/
	mv ${wwwdir}/topsel/unweighted/*.{png,pdf} ${wwwdir}/topsel/unweighted/norm/
	mv ${wwwdir}/topsel/unweighted/norm/*_all*.{png,pdf} ${wwwdir}/topsel/unweighted/norm/all/
	mv ${wwwdir}/topsel/unweighted/norm/*_ee*.{png,pdf} ${wwwdir}/topsel/unweighted/norm/ee/
	mv ${wwwdir}/topsel/unweighted/norm/*_em*.{png,pdf} ${wwwdir}/topsel/unweighted/norm/em/
	mv ${wwwdir}/topsel/unweighted/norm/*_mm*.{png,pdf} ${wwwdir}/topsel/unweighted/norm/mm/
	mv ${wwwdir}/topsel/unweighted/norm/*_e*.{png,pdf} ${wwwdir}/topsel/unweighted/norm/e/
	mv ${wwwdir}/topsel/unweighted/norm/*_m*.{png,pdf} ${wwwdir}/topsel/unweighted/norm/m/ 
  mv ${wwwdir}/topsel/*.{png,pdf} ${wwwdir}/topsel/weighted/
	mv ${wwwdir}/topsel/weighted/*_log*.{png,pdf} ${wwwdir}/topsel/weighted/log/
	mv ${wwwdir}/topsel/weighted/log/*_all*.{png,pdf} ${wwwdir}/topsel/weighted/log/all/
	mv ${wwwdir}/topsel/weighted/log/*_ee*.{png,pdf} ${wwwdir}/topsel/weighted/log/ee/
	mv ${wwwdir}/topsel/weighted/log/*_em*.{png,pdf} ${wwwdir}/topsel/weighted/log/em/
	mv ${wwwdir}/topsel/weighted/log/*_mm*.{png,pdf} ${wwwdir}/topsel/weighted/log/mm/
	mv ${wwwdir}/topsel/weighted/log/*_e*.{png,pdf} ${wwwdir}/topsel/weighted/log/e/
	mv ${wwwdir}/topsel/weighted/log/*_m*.{png,pdf} ${wwwdir}/topsel/weighted/log/m/
	mv ${wwwdir}/topsel/weighted/*.{png,pdf} ${wwwdir}/topsel/weighted/norm/
	mv ${wwwdir}/topsel/weighted/norm/*_all*.{png,pdf} ${wwwdir}/topsel/weighted/norm/all/
	mv ${wwwdir}/topsel/weighted/norm/*_ee*.{png,pdf} ${wwwdir}/topsel/weighted/norm/ee/
	mv ${wwwdir}/topsel/weighted/norm/*_em*.{png,pdf} ${wwwdir}/topsel/weighted/norm/em/
	mv ${wwwdir}/topsel/weighted/norm/*_mm*.{png,pdf} ${wwwdir}/topsel/weighted/norm/mm/
	mv ${wwwdir}/topsel/weighted/norm/*_e*.{png,pdf} ${wwwdir}/topsel/weighted/norm/e/
	mv ${wwwdir}/topsel/weighted/norm/*_m*.{png,pdf} ${wwwdir}/topsel/weighted/norm/m/ 
   
	mv ${wwwdir}/*_jpsicand*.{png,pdf} ${wwwdir}/jpsicand/
	mv ${wwwdir}/jpsicand/*_no_weight*.{png,pdf} ${wwwdir}/jpsicand/unweighted/
	mv ${wwwdir}/jpsicand/unweighted/*_log*.{png,pdf} ${wwwdir}/jpsicand/unweighted/log/
	mv ${wwwdir}/jpsicand/unweighted/log/*_all*.{png,pdf} ${wwwdir}/jpsicand/unweighted/log/all/
	mv ${wwwdir}/jpsicand/unweighted/log/*_ee*.{png,pdf} ${wwwdir}/jpsicand/unweighted/log/ee/
	mv ${wwwdir}/jpsicand/unweighted/log/*_em*.{png,pdf} ${wwwdir}/jpsicand/unweighted/log/em/
	mv ${wwwdir}/jpsicand/unweighted/log/*_mm*.{png,pdf} ${wwwdir}/jpsicand/unweighted/log/mm/
	mv ${wwwdir}/jpsicand/unweighted/log/*_e*.{png,pdf} ${wwwdir}/jpsicand/unweighted/log/e/
	mv ${wwwdir}/jpsicand/unweighted/log/*_m*.{png,pdf} ${wwwdir}/jpsicand/unweighted/log/m/
	mv ${wwwdir}/jpsicand/unweighted/*.{png,pdf} ${wwwdir}/jpsicand/unweighted/norm/
	mv ${wwwdir}/jpsicand/unweighted/norm/*_all*.{png,pdf} ${wwwdir}/jpsicand/unweighted/norm/all/
	mv ${wwwdir}/jpsicand/unweighted/norm/*_ee*.{png,pdf} ${wwwdir}/jpsicand/unweighted/norm/ee/
	mv ${wwwdir}/jpsicand/unweighted/norm/*_em*.{png,pdf} ${wwwdir}/jpsicand/unweighted/norm/em/
	mv ${wwwdir}/jpsicand/unweighted/norm/*_mm*.{png,pdf} ${wwwdir}/jpsicand/unweighted/norm/mm/
	mv ${wwwdir}/jpsicand/unweighted/norm/*_e*.{png,pdf} ${wwwdir}/jpsicand/unweighted/norm/e/
	mv ${wwwdir}/jpsicand/unweighted/norm/*_m*.{png,pdf} ${wwwdir}/jpsicand/unweighted/norm/m/ 
  mv ${wwwdir}/jpsicand/*.{png,pdf} ${wwwdir}/jpsicand/weighted/
	mv ${wwwdir}/jpsicand/weighted/*_log*.{png,pdf} ${wwwdir}/jpsicand/weighted/log/
	mv ${wwwdir}/jpsicand/weighted/log/*_all*.{png,pdf} ${wwwdir}/jpsicand/weighted/log/all/
	mv ${wwwdir}/jpsicand/weighted/log/*_ee*.{png,pdf} ${wwwdir}/jpsicand/weighted/log/ee/
	mv ${wwwdir}/jpsicand/weighted/log/*_em*.{png,pdf} ${wwwdir}/jpsicand/weighted/log/em/
	mv ${wwwdir}/jpsicand/weighted/log/*_mm*.{png,pdf} ${wwwdir}/jpsicand/weighted/log/mm/
	mv ${wwwdir}/jpsicand/weighted/log/*_e*.{png,pdf} ${wwwdir}/jpsicand/weighted/log/e/
	mv ${wwwdir}/jpsicand/weighted/log/*_m*.{png,pdf} ${wwwdir}/jpsicand/weighted/log/m/
	mv ${wwwdir}/jpsicand/weighted/*.{png,pdf} ${wwwdir}/jpsicand/weighted/norm/
	mv ${wwwdir}/jpsicand/weighted/norm/*_all*.{png,pdf} ${wwwdir}/jpsicand/weighted/norm/all/
	mv ${wwwdir}/jpsicand/weighted/norm/*_ee*.{png,pdf} ${wwwdir}/jpsicand/weighted/norm/ee/
	mv ${wwwdir}/jpsicand/weighted/norm/*_em*.{png,pdf} ${wwwdir}/jpsicand/weighted/norm/em/
	mv ${wwwdir}/jpsicand/weighted/norm/*_mm*.{png,pdf} ${wwwdir}/jpsicand/weighted/norm/mm/
	mv ${wwwdir}/jpsicand/weighted/norm/*_e*.{png,pdf} ${wwwdir}/jpsicand/weighted/norm/e/
	mv ${wwwdir}/jpsicand/weighted/norm/*_m*.{png,pdf} ${wwwdir}/jpsicand/weighted/norm/m/ 
	
	mv ${wwwdir}/*_d0cand*.{png,pdf} ${wwwdir}/d0cand/
	mv ${wwwdir}/d0cand/*_no_weight*.{png,pdf} ${wwwdir}/d0cand/unweighted/
	mv ${wwwdir}/d0cand/unweighted/*_log*.{png,pdf} ${wwwdir}/d0cand/unweighted/log/
	mv ${wwwdir}/d0cand/unweighted/log/*_all*.{png,pdf} ${wwwdir}/d0cand/unweighted/log/all/
	mv ${wwwdir}/d0cand/unweighted/log/*_ee*.{png,pdf} ${wwwdir}/d0cand/unweighted/log/ee/
	mv ${wwwdir}/d0cand/unweighted/log/*_em*.{png,pdf} ${wwwdir}/d0cand/unweighted/log/em/
	mv ${wwwdir}/d0cand/unweighted/log/*_mm*.{png,pdf} ${wwwdir}/d0cand/unweighted/log/mm/
	mv ${wwwdir}/d0cand/unweighted/log/*_e*.{png,pdf} ${wwwdir}/d0cand/unweighted/log/e/
	mv ${wwwdir}/d0cand/unweighted/log/*_m*.{png,pdf} ${wwwdir}/d0cand/unweighted/log/m/
	mv ${wwwdir}/d0cand/unweighted/*.{png,pdf} ${wwwdir}/d0cand/unweighted/norm/
	mv ${wwwdir}/d0cand/unweighted/norm/*_all*.{png,pdf} ${wwwdir}/d0cand/unweighted/norm/all/
	mv ${wwwdir}/d0cand/unweighted/norm/*_ee*.{png,pdf} ${wwwdir}/d0cand/unweighted/norm/ee/
	mv ${wwwdir}/d0cand/unweighted/norm/*_em*.{png,pdf} ${wwwdir}/d0cand/unweighted/norm/em/
	mv ${wwwdir}/d0cand/unweighted/norm/*_mm*.{png,pdf} ${wwwdir}/d0cand/unweighted/norm/mm/
	mv ${wwwdir}/d0cand/unweighted/norm/*_e*.{png,pdf} ${wwwdir}/d0cand/unweighted/norm/e/
	mv ${wwwdir}/d0cand/unweighted/norm/*_m*.{png,pdf} ${wwwdir}/d0cand/unweighted/norm/m/ 
  mv ${wwwdir}/d0cand/*.{png,pdf} ${wwwdir}/d0cand/weighted/
	mv ${wwwdir}/d0cand/weighted/*_log*.{png,pdf} ${wwwdir}/d0cand/weighted/log/
	mv ${wwwdir}/d0cand/weighted/log/*_all*.{png,pdf} ${wwwdir}/d0cand/weighted/log/all/
	mv ${wwwdir}/d0cand/weighted/log/*_ee*.{png,pdf} ${wwwdir}/d0cand/weighted/log/ee/
	mv ${wwwdir}/d0cand/weighted/log/*_em*.{png,pdf} ${wwwdir}/d0cand/weighted/log/em/
	mv ${wwwdir}/d0cand/weighted/log/*_mm*.{png,pdf} ${wwwdir}/d0cand/weighted/log/mm/
	mv ${wwwdir}/d0cand/weighted/log/*_e*.{png,pdf} ${wwwdir}/d0cand/weighted/log/e/
	mv ${wwwdir}/d0cand/weighted/log/*_m*.{png,pdf} ${wwwdir}/d0cand/weighted/log/m/
	mv ${wwwdir}/d0cand/weighted/*.{png,pdf} ${wwwdir}/d0cand/weighted/norm/
	mv ${wwwdir}/d0cand/weighted/norm/*_all*.{png,pdf} ${wwwdir}/d0cand/weighted/norm/all/
	mv ${wwwdir}/d0cand/weighted/norm/*_ee*.{png,pdf} ${wwwdir}/d0cand/weighted/norm/ee/
	mv ${wwwdir}/d0cand/weighted/norm/*_em*.{png,pdf} ${wwwdir}/d0cand/weighted/norm/em/
	mv ${wwwdir}/d0cand/weighted/norm/*_mm*.{png,pdf} ${wwwdir}/d0cand/weighted/norm/mm/
	mv ${wwwdir}/d0cand/weighted/norm/*_e*.{png,pdf} ${wwwdir}/d0cand/weighted/norm/e/
	mv ${wwwdir}/d0cand/weighted/norm/*_m*.{png,pdf} ${wwwdir}/d0cand/weighted/norm/m/ 

	cp -p test/index.php ${wwwdir}/
	
  cp -p test/index.php ${wwwdir}/check/
	cp -p test/index.php ${wwwdir}/check/unweighted/
	cp -p test/index.php ${wwwdir}/check/unweighted/log/
	cp -p test/index.php ${wwwdir}/check/unweighted/log/all/
	cp -p test/index.php ${wwwdir}/check/unweighted/log/ee/
	cp -p test/index.php ${wwwdir}/check/unweighted/log/em/
	cp -p test/index.php ${wwwdir}/check/unweighted/log/mm/
	cp -p test/index.php ${wwwdir}/check/unweighted/log/e/
	cp -p test/index.php ${wwwdir}/check/unweighted/log/m/
	cp -p test/index.php ${wwwdir}/check/unweighted/norm/
	cp -p test/index.php ${wwwdir}/check/unweighted/norm/all/
	cp -p test/index.php ${wwwdir}/check/unweighted/norm/ee/
	cp -p test/index.php ${wwwdir}/check/unweighted/norm/em/
	cp -p test/index.php ${wwwdir}/check/unweighted/norm/mm/
	cp -p test/index.php ${wwwdir}/check/unweighted/norm/e/
	cp -p test/index.php ${wwwdir}/check/unweighted/norm/m/
	cp -p test/index.php ${wwwdir}/check/weighted/
	cp -p test/index.php ${wwwdir}/check/weighted/log/
	cp -p test/index.php ${wwwdir}/check/weighted/log/all/
	cp -p test/index.php ${wwwdir}/check/weighted/log/ee/
	cp -p test/index.php ${wwwdir}/check/weighted/log/em/
	cp -p test/index.php ${wwwdir}/check/weighted/log/mm/
	cp -p test/index.php ${wwwdir}/check/weighted/log/e/
	cp -p test/index.php ${wwwdir}/check/weighted/log/m/
	cp -p test/index.php ${wwwdir}/check/weighted/norm/
	cp -p test/index.php ${wwwdir}/check/weighted/norm/all/
	cp -p test/index.php ${wwwdir}/check/weighted/norm/ee/
	cp -p test/index.php ${wwwdir}/check/weighted/norm/em/
	cp -p test/index.php ${wwwdir}/check/weighted/norm/mm/
	cp -p test/index.php ${wwwdir}/check/weighted/norm/e/
	cp -p test/index.php ${wwwdir}/check/weighted/norm/m/
	
  cp -p test/index.php ${wwwdir}/topsel/
	cp -p test/index.php ${wwwdir}/topsel/unweighted/
	cp -p test/index.php ${wwwdir}/topsel/unweighted/log/
	cp -p test/index.php ${wwwdir}/topsel/unweighted/log/all/
	cp -p test/index.php ${wwwdir}/topsel/unweighted/log/ee/
	cp -p test/index.php ${wwwdir}/topsel/unweighted/log/em/
	cp -p test/index.php ${wwwdir}/topsel/unweighted/log/mm/
	cp -p test/index.php ${wwwdir}/topsel/unweighted/log/e/
	cp -p test/index.php ${wwwdir}/topsel/unweighted/log/m/
	cp -p test/index.php ${wwwdir}/topsel/unweighted/norm/
	cp -p test/index.php ${wwwdir}/topsel/unweighted/norm/all/
	cp -p test/index.php ${wwwdir}/topsel/unweighted/norm/ee/
	cp -p test/index.php ${wwwdir}/topsel/unweighted/norm/em/
	cp -p test/index.php ${wwwdir}/topsel/unweighted/norm/mm/
	cp -p test/index.php ${wwwdir}/topsel/unweighted/norm/e/
	cp -p test/index.php ${wwwdir}/topsel/unweighted/norm/m/
	cp -p test/index.php ${wwwdir}/topsel/weighted/
	cp -p test/index.php ${wwwdir}/topsel/weighted/log/
	cp -p test/index.php ${wwwdir}/topsel/weighted/log/all/
	cp -p test/index.php ${wwwdir}/topsel/weighted/log/ee/
	cp -p test/index.php ${wwwdir}/topsel/weighted/log/em/
	cp -p test/index.php ${wwwdir}/topsel/weighted/log/mm/
	cp -p test/index.php ${wwwdir}/topsel/weighted/log/e/
	cp -p test/index.php ${wwwdir}/topsel/weighted/log/m/
	cp -p test/index.php ${wwwdir}/topsel/weighted/norm/
	cp -p test/index.php ${wwwdir}/topsel/weighted/norm/all/
	cp -p test/index.php ${wwwdir}/topsel/weighted/norm/ee/
	cp -p test/index.php ${wwwdir}/topsel/weighted/norm/em/
	cp -p test/index.php ${wwwdir}/topsel/weighted/norm/mm/
	cp -p test/index.php ${wwwdir}/topsel/weighted/norm/e/
	cp -p test/index.php ${wwwdir}/topsel/weighted/norm/m/
	
  cp -p test/index.php ${wwwdir}/jpsicand/
	cp -p test/index.php ${wwwdir}/jpsicand/unweighted/
	cp -p test/index.php ${wwwdir}/jpsicand/unweighted/log/
	cp -p test/index.php ${wwwdir}/jpsicand/unweighted/log/all/
	cp -p test/index.php ${wwwdir}/jpsicand/unweighted/log/ee/
	cp -p test/index.php ${wwwdir}/jpsicand/unweighted/log/em/
	cp -p test/index.php ${wwwdir}/jpsicand/unweighted/log/mm/
	cp -p test/index.php ${wwwdir}/jpsicand/unweighted/log/e/
	cp -p test/index.php ${wwwdir}/jpsicand/unweighted/log/m/
	cp -p test/index.php ${wwwdir}/jpsicand/unweighted/norm/
	cp -p test/index.php ${wwwdir}/jpsicand/unweighted/norm/all/
	cp -p test/index.php ${wwwdir}/jpsicand/unweighted/norm/ee/
	cp -p test/index.php ${wwwdir}/jpsicand/unweighted/norm/em/
	cp -p test/index.php ${wwwdir}/jpsicand/unweighted/norm/mm/
	cp -p test/index.php ${wwwdir}/jpsicand/unweighted/norm/e/
	cp -p test/index.php ${wwwdir}/jpsicand/unweighted/norm/m/
	cp -p test/index.php ${wwwdir}/jpsicand/weighted/
	cp -p test/index.php ${wwwdir}/jpsicand/weighted/log/
	cp -p test/index.php ${wwwdir}/jpsicand/weighted/log/all/
	cp -p test/index.php ${wwwdir}/jpsicand/weighted/log/ee/
	cp -p test/index.php ${wwwdir}/jpsicand/weighted/log/em/
	cp -p test/index.php ${wwwdir}/jpsicand/weighted/log/mm/
	cp -p test/index.php ${wwwdir}/jpsicand/weighted/log/e/
	cp -p test/index.php ${wwwdir}/jpsicand/weighted/log/m/
	cp -p test/index.php ${wwwdir}/jpsicand/weighted/norm/
	cp -p test/index.php ${wwwdir}/jpsicand/weighted/norm/all/
	cp -p test/index.php ${wwwdir}/jpsicand/weighted/norm/ee/
	cp -p test/index.php ${wwwdir}/jpsicand/weighted/norm/em/
	cp -p test/index.php ${wwwdir}/jpsicand/weighted/norm/mm/
	cp -p test/index.php ${wwwdir}/jpsicand/weighted/norm/e/
	cp -p test/index.php ${wwwdir}/jpsicand/weighted/norm/m/

	cp -p test/index.php ${wwwdir}/d0cand/
	cp -p test/index.php ${wwwdir}/d0cand/unweighted/
	cp -p test/index.php ${wwwdir}/d0cand/unweighted/log/
	cp -p test/index.php ${wwwdir}/d0cand/unweighted/log/all/
	cp -p test/index.php ${wwwdir}/d0cand/unweighted/log/ee/
	cp -p test/index.php ${wwwdir}/d0cand/unweighted/log/em/
	cp -p test/index.php ${wwwdir}/d0cand/unweighted/log/mm/
	cp -p test/index.php ${wwwdir}/d0cand/unweighted/log/e/
	cp -p test/index.php ${wwwdir}/d0cand/unweighted/log/m/
	cp -p test/index.php ${wwwdir}/d0cand/unweighted/norm/
	cp -p test/index.php ${wwwdir}/d0cand/unweighted/norm/all/
	cp -p test/index.php ${wwwdir}/d0cand/unweighted/norm/ee/
	cp -p test/index.php ${wwwdir}/d0cand/unweighted/norm/em/
	cp -p test/index.php ${wwwdir}/d0cand/unweighted/norm/mm/
	cp -p test/index.php ${wwwdir}/d0cand/unweighted/norm/e/
	cp -p test/index.php ${wwwdir}/d0cand/unweighted/norm/m/
	cp -p test/index.php ${wwwdir}/d0cand/weighted/
	cp -p test/index.php ${wwwdir}/d0cand/weighted/log/
	cp -p test/index.php ${wwwdir}/d0cand/weighted/log/all/
	cp -p test/index.php ${wwwdir}/d0cand/weighted/log/ee/
	cp -p test/index.php ${wwwdir}/d0cand/weighted/log/em/
	cp -p test/index.php ${wwwdir}/d0cand/weighted/log/mm/
	cp -p test/index.php ${wwwdir}/d0cand/weighted/log/e/
	cp -p test/index.php ${wwwdir}/d0cand/weighted/log/m/
	cp -p test/index.php ${wwwdir}/d0cand/weighted/norm/
	cp -p test/index.php ${wwwdir}/d0cand/weighted/norm/all/
	cp -p test/index.php ${wwwdir}/d0cand/weighted/norm/ee/
	cp -p test/index.php ${wwwdir}/d0cand/weighted/norm/em/
	cp -p test/index.php ${wwwdir}/d0cand/weighted/norm/mm/
	cp -p test/index.php ${wwwdir}/d0cand/weighted/norm/e/
	cp -p test/index.php ${wwwdir}/d0cand/weighted/norm/m/

  echo "Please visit http://ebouvier.web.cern.ch/ebouvier/Elvire/"
	;;
esac
