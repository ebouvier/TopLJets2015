#!/bin/bash

WHAT=$1; 
if [[ "$1" == "" ]]; then 
    echo "steerTOPWidthAnalysis.sh <SEL/MERGESEL/PLOTSEL/WWWSEL/ANA/MERGE/PLOT/WWW>";
    echo "        SEL        - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGESEL   - merge the output of the jobs";
    echo "        PLOTSEL    - runs the plotter tool on the selection";
    echo "        WWWSEL     - moves the plots to an afs-web-based area";
    echo "        ANA        - analyze the selected events";
    echo "        MERGE      - merge the output of the analysis jobs";
    echo "        PLOT       - runs the plotter tool on the analysis outputs";
    echo "        WWW        - moves the analysis plots to an afs-web-based area";

    exit 1; 
fi

queue=8nh
eosdir=/store/cmst3/user/psilva/LJets2015/8c1e7c9
#outdir=~/work/TopWidth
outdir=/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/analysis
wwwdir=~/www/Top2016
method=TOP::RunTop
lumi=2267.84

RED='\e[31m'
NC='\e[0m'

case $WHAT in
    SEL )
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${outdir} -m ${method} --ch 0;
	;;
    MERGESEL )
	./scripts/mergeOutputs.py ${outdir};	
	;;
    PLOTSEL )
	#python scripts/plotter.py -i ${outdir} --puNormSF puwgtctr  -j data/samples_Run2015.json -l ${lumi};	
	python scripts/plotter.py -i ${outdir} -j data/samples_Run2015.json -l ${lumi};	
	;;
    WWWSEL )
	mkdir -p ${wwwdir}/sel
	cp ${outdir}/plots/*.{png,pdf} ${wwwdir}/sel
	cp test/index.php ${wwwdir}/sel
	;;
    ANA )
	python scripts/runTopWidthAnalysis.py -i ${outdir}/Chunks -o ${outdir}/analysis -q 8nh;
	;;
    MERGE )
	./scripts/mergeOutputs.py ${outdir}/analysis;
	;;
    PLOT )
        python scripts/plotter.py -i ${outdir}/analysis  -j data/samples_Run2015.json -l ${lumi};        
        ;;
    WWW )
        mkdir -p ${wwwdir}/ana
        cp ${outdir}/plots/*.{png,pdf} ${wwwdir}/ana        
        cp test/index.php ${wwwdir}/ana
	;;
esac
