#!/bin/bash
queue=8nh
eosdir=/store/cmst3/user/psilva/LJets2015/8c1e7c9
#outdir=/store/user/byates/LJets2015/skim
outdir=/afs/cern.ch/user/b/byates/eos/cms/store/user/byates/LJets2015/skim
wwwdir=~/www/LJets2015-arcrev
lumi=2267.84
#curr_dir=`pwd`
#ANAL_DIR=/afs/cern.ch/user/b/byates/TopAnalysis
#cd $ANAL_DIR
#eval `scram runtime -sh`
#cd $curr_dir

#csv=(Data13TeV_SingleElectron_2015C,Data13TeV_SingleElectron_2015D,Data13TeV_SingleMuon_2015C,Data13TeV_SingleMuon_2015D)

#python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} --runSysts -o ${outdir}/analysis_muplus  --ch 13   --charge 1 --only ${csv}
#python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} --runSysts -o ${outdir}/analysis_muminus  --ch 13   --charge -1 --only ${csv}
#python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} --runSysts -o ${outdir}/analysis_eplus   --ch 11   --charge 1 --only ${csv}
#python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} --runSysts -o ${outdir}/analysis_eminus  --ch 11   --charge -1 --only ${csv}

python scripts/runLocalAnalysis.py -i ${eosdir} -o ${outdir} -q ${queue} -n 8 --method TOPWidth::RunTopWidth
#mv *.root $outdir/
