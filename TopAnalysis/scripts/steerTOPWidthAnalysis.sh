#!/bin/bash

WHAT=$1; 
ERA=$2
if [ "$#" -ne 2 ]; then 
    echo "steerTOPWidthAnalysis.sh <SEL/MERGESEL/PLOTSEL/WWWSEL/ANA/MERGE/BKG/PLOT/WWW> <ERA>";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGESEL     - merge the output of the jobs";
    echo "        PLOTSEL      - runs the plotter tool on the selection";
    echo "        WWWSEL       - moves the plots to an afs-web-based area";
    echo "        ANA          - analyze the selected events";
    echo "        MERGE        - merge the output of the analysis jobs";
    echo "        BKG          - estimate DY scale factor from data";
    echo "        PLOT         - runs the plotter tool on the analysis outputs";
    echo "        WWW          - moves the analysis plots to an afs-web-based area";
    echo " "
    echo "        ERA          - era2015/era2016";
    exit 1; 
fi

export LSB_JOB_REPORT_MAIL=N

queue=8nh
githash=8db9ad6
#lumi=3977.28
lumi=12868.66
lumiUnc=0.062
#eosdir=/store/cmst3/user/psilva/LJets2016/${githash}
eosdir=/store/user/byates/LJets2015/${githash}
case $ERA in
    era2015)
	githash=8c1e7c9;
	lumi=2267.84
	eosdir=/store/cmst3/user/psilva/LJets2015/${githash}
	;;
esac

summaryeosdir=/store/cmst3/group/top/summer2016/TopWidth_${ERA}
#outdir=~/work/TopWidth_${ERA}
outdir=/afs/cern.ch/user/b/byates/CMSSW_8_0_11/src/TopLJets2015/TopAnalysis/LJets2015/2016
wwwdir=~/www/Top2016/2016


RED='\e[31m'
NC='\e[0m'
case $WHAT in
    SEL )
	python scripts/runLocalAnalysis.py -i ${eosdir} -n 8 -q ${queue} -o ${outdir} --era ${ERA} -m TOP::RunTop --ch 0;
	;;
    MERGESEL )
	./scripts/mergeOutputs.py ${outdir} True;	
	;;
    PLOTSEL )
	python scripts/plotter.py -i ${outdir} --puNormSF puwgtctr  -j data/${ERA}/samples.json -l ${lumi} --saveLog;# --mcUnc ${lumiUnc};	
        #python scripts/plotter.py -i ${outdir} -j data/${ERA}/samples.json -l ${lumi};
	;;
    WWWSEL )
	mkdir -p ${wwwdir}/sel
	cp -p  ${outdir}/plots/*.{png,pdf} ${wwwdir}/sel/$subdir
	rm ${wwwdir}/sel/*no_weight*.{png,pdf}
	mv ${wwwdir}/sel/*_log.{png,pdf} ${wwwdir}/sel/log/
	mv ${wwwdir}/sel/*_ee*.{png,pdf} ${wwwdir}/sel/ee/
	mv ${wwwdir}/sel/*_em*.{png,pdf} ${wwwdir}/sel/em/
	mv ${wwwdir}/sel/*_e*.{png,pdf} ${wwwdir}/sel/e/
	mv ${wwwdir}/sel/*_mm*.{png,pdf} ${wwwdir}/sel/mumu/
	mv ${wwwdir}/sel/*_m*.{png,pdf} ${wwwdir}/sel/mu/
	cp -p test/index.php ${wwwdir}/sel/$subdir
	cp -p test/index.php ${wwwdir}/sel/log/
	cp -p test/index.php ${wwwdir}/sel/ee/
	cp -p test/index.php ${wwwdir}/sel/e/
	cp -p test/index.php ${wwwdir}/sel/em/
	cp -p test/index.php ${wwwdir}/sel/mumu/
	cp -p test/index.php ${wwwdir}/sel/mu/
	;;
    ANA )
	python scripts/runTopWidthAnalysis.py -i ${summaryeosdir} -o ${outdir}/analysis -q ${queue};
	;;
    MERGE )
	./scripts/mergeOutputs.py ${outdir}/analysis;
	;;
    BKG )
	python scripts/plotter.py      -i ${outdir}/analysis  -j data/${ERA}/samples.json  -l ${lumi} --onlyData --only mll -o dy_plotter.root;        
	python scripts/runDYRinRout.py --in ${outdir}/analysis/plots/dy_plotter.root --categs 1b,2b --out ${outdir}/analysis/plots/;
	;;
    PLOT )
	python scripts/plotter.py -i ${outdir}/analysis  -j data/${ERA}/samples.json      -l ${lumi} --only count --saveTeX -o count_plotter.root --procSF DY:${outdir}/analysis/plots/.dyscalefactors.pck; 
        python scripts/plotter.py -i ${outdir}/analysis  -j data/${ERA}/samples.json      -l ${lumi} --onlyData --procSF DY:${outdir}/analysis/plots/.dyscalefactors.pck;        
	python scripts/plotter.py -i ${outdir}/analysis  -j data/${ERA}/syst_samples.json -l ${lumi} --silent -o syst_plotter.root;        
        ;;
    WWW )
        mkdir -p ${wwwdir}/ana
        cp ${outdir}/analysis/plots/*.{png,pdf} ${wwwdir}/ana        
        cp test/index.php ${wwwdir}/ana
	;;
esac
