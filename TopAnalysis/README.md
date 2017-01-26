# TopLJets2015

## References

[Official TWiki](https://twiki.cern.ch/twiki/bin/view/Main/TopLJ2015Analysis), 
[upstream Git repository](https://github.com/bryates/TopLJets2015.git), and
[associated plots](https://byates.web.cern.ch/byates/Top2016/2016/sel/).

## Installation instructions

To execute in your lxplus work area.
```
cmsrel CMSSW_8_0_11 
cd CMSSW_8_0_11/src 
cmsenv
git cms-init
#EGM smearer
git remote add -f -t ecal_smear_fix_80X emanueledimarco https://github.com/emanueledimarco/cmssw.git
git cms-addpkg EgammaAnalysis/ElectronTools
git checkout -b from-52f192a 52f192a
cd EgammaAnalysis/ElectronTools/data
git clone -b ICHEP2016_v2 https://github.com/ECALELFS/ScalesSmearings.git
cd -
#Pseudo-top producer with Markus fix
git cms-addpkg  TopQuarkAnalysis/TopEventProducers
https://raw.githubusercontent.com/intrepid42/cmssw/4336e8182cab054c8383d7b4eb6622c046952711/TopQuarkAnalysis/TopEventProducers/src/PseudoTopProducer.cc
#analysis code
cd-
git clone git@github.com:ebouvier/TopLJets2015.git
cd TopLJets2015/TopAnalysis
git checkout 80x_dev
scram b -j 8
```

## Running ntuple creation

To run locally the ntuplizer, for testing purposes
```
cmsRun test/runMiniAnalyzer_cfg.py runOnData=False/True outFilename=MiniEvents.root
```
To submit a list of samples, described in a json file to the grid you can use the following script.
```
python scripts/submitToGrid.py -j data/era2016/samples.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py --lfn /store/group/phys_top/byates -s
```
Partial submission can be made adding "-o csv_list" as an option
Don't forget to init the environment for crab3 (e.g. https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3CheatSheet#Environment_setup)
```
source /cvmfs/cms.cern.ch/crab3/crab.sh
```
As soon as ntuple production starts to finish, to move from crab output directories to a simpler directory structure which can be easily parsed by the local analysis runThe merging can be run locally if needed by using the checkProductionIntegrity.py script
```
python scripts/submitCheckProductionIntegrity.py -i /store/group/phys_top/byates/65ee28f -o /store/user/byates/LJets2015/8db9ad6
```

## Preparing the analysis 

Correction and uncertainty files are stored under data by era directories (e.g. data/era2015, data/era2016) in order no to mix different periods.
After ntuples are processed start by creating the json files with the list of runs/luminosity sections processed, e.g. as:
```
crab report grid/crab_Data13TeV_DoubleMuon_2016B
``` 
Then you can merge the json files for the same dataset to get the full list of run/lumi sections to analyse
```
mergeJSON.py grid/crab_Data13TeV_DoubleMuon_2016B/results/processedLumis.json grid/crab_Data13TeV_DoubleMuon_2015C/results/processedLumis.json grid/crab_Data13TeV_DoubleMuon_2015D/results/processedLumis.json --output data/era2016/Data13TeV_DoubleMuon_lumis.json
```
You can then run the brilcalc tool to get the integrated luminosity in total and per run (see https://twiki.cern.ch/twiki/bin/view/CMS/2015LumiNormtag for more details).
```
export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.0.3/bin:$PATH
brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i data/era2016/Data13TeV_DoubleMuon_lumis.json
```
Use the table which is printed out to update the "lumiPerRun" method in ReadTree.cc.
That will be used to monitor the event yields per run in order to identify outlier runs.
* Pileup weighting. To update the pileup distributions run the script below. It will store the data pileup distributions for different min.bias cross section in data/pileupWgts.root
```
python scripts/runPileupEstimation.py --json data/era2016/Data13TeV_DoubleMuon_lumis.json --out data/era2016/pileupWgts.root
```
* B-tagging. To apply corrections to the simulation one needs the expected efficiencies stored somwewhere. The script below will project the jet pT spectrum from the TTbar sample before and after applying b-tagging, to compute the expecte efficiencies. The result will be stored in data/expTageff.root
```
for i in "_powheg" "_herwig" "_scaledown" "_scaleup"; do
    python scripts/saveExpectedBtagEff.py -i /store/user/byates/LJets2015/8db9ad6/MC13TeV_TTJets${i} -o data/era2016/expTageff${i}.root;
done
mv data/era2016/expTageff_powheg.root data/era2016/expTageff.root
```
* MC normalization. This will loop over all the samples available in EOS and produce a normalization cache (weights to normalize MC). The file will be available in data/genweights.pck
```
python scripts/produceNormalizationCache.py -i /store/user/byates/LJets2015/8db9ad6 -o data/era2016/genweights.root
```
You're now ready to start locally the analysis.


## Running locally the analysis for testing

The analysis (histogram filling, final selection) is in src/ReadTree.cc.
Recompile (scram b) everytime you change it so that you can test the new features.
To test the code on a single file to produce plots.
```
python scripts/runLocalAnalysis.py -i MiniEvents.root
```
To run the code on a set of samples stored in EOS you can run it as shown below.
If "-q queue_name" is appended the jobs are submitted to the batch system instead of running locally. 
To check the status of your jobs run "bjobs" and then "bpeek job_number" if you want to inspect how the job is running in the cluster.
If "-n n_jobs" is passed the script runs locally using "n_jobs" parallel threads.
```
python scripts/runLocalAnalysis.py -i /store/user/byates/LJets2015/8db9ad6 -n 8 --runSysts -o analysis_muplus   --ch 13   --charge 1
```
If you want to suppress the mails sent automatically after job completion please do
```
export LSB_JOB_REPORT_MAIL=N
```
before submitting the jobs to the batch. After the jobs have run you can merge the outputs with
```
./scripts/mergeOutputs.py analysis_muplus
```
To plot the output of the local analysis you can run the following:
```
python scripts/plotter.py -i analysis_muplus/   -j data/era2016/samples.json  -l 12870
```

## Submitting the full analysis to the batch system

A script wraps up the above procedure for all the signal and control regions used in the analyis.
To use it you can use the following script
```
sh scripts/steerTOPMassAnalysis.sh <SEL/MERGESEL/PLOTSEL/WWWSEL>
```
Job status can be check with the command `bjobs`. For more details about the batch system, see the [dedicated TWiki page](https://twiki.cern.ch/twiki/bin/view/Main/BatchJobs).
Finally, to get the detailed event yields, run:
```
python getNumberOfEvents.py
```


## Updating the code

Commit your changes regularly with
```
git commit -a -m'comment on the changes made'
```
Push to your forked repository
```
git push git@github.com:MYGITHUBLOGIN/TopLJets2015.git
```
From the github area of the repository cleak on the green button "Compare,review and create a pull request"
to create the PR to merge with your colleagues.
