#!/bin/bash
queue=8nh
eosdir=/store/cmst3/user/psilva/LJets2015/8c1e7c9
outdir=/afs/cern.ch/user/b/byates/TopAnalysis/LJets2015/skim
wwwdir=~/www/LJets2015-arcrev
lumi=2267.84
curr_dir=`pwd`
#ANAL_DIR=/afs/cern.ch/user/b/byates/TopAnalysis
#cd $ANAL_DIR
#eval `scram runtime -sh`
#cd $curr_dir

#csv=MC13TeV_WJets_madgraph,MC13TeV_WWTo2L2Nu,MC13TeV_WWToLNuQQ,MC13TeV_WZ,MC13TeV_ZZ,MC13TeV_tW_DS,MC13TeV_tW_m169v5,MC13TeV_tW_m175v5,MC13TeV_tW_scaledown,MC13TeV_tW_scaleup,MC13TeV_tbarW_DS,MC13TeV_tbarW_m169v5,MC13TeV_tbarW_m175v5,MC13TeV_tbarW_scaledown,MC13TeV_tbarW_scaleup
#csv=MC13TeV_TTJets_MPIoff,MC13TeV_TTJets_herwig,MC13TeV_TTJets_m169v5,MC13TeV_TTJets_m175v5,MC13TeV_TTJets_noCR,MC13TeV_TTJets_scaledown,MC13TeV_TTJets_scaleup,MC13TeV_TTWToLNu,MC13TeV_TTWToQQ,MC13TeV_TTZToLLNuNu,MC13TeV_TTZToQQ
csv=MC13TeV_TTJets_amcatnloFXFX
#python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} --runSysts -o ${outdir}/analysis_muplus  --ch 13   --charge 1 --only ${csv}
#python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} --runSysts -o ${outdir}/analysis_muminus  --ch 13   --charge -1 --only ${csv}
#python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} --runSysts -o ${outdir}/analysis_eplus   --ch 11   --charge 1 --only ${csv}
#python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} --runSysts -o ${outdir}/analysis_eminus  --ch 11   --charge -1 --only ${csv}

python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${outdir} --method TOPWidth::RunTopWidth --only ${csv}
