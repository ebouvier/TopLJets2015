#!/bin/bash
for i in "" "_herwig" "_scaledown" "_scaleup"; do
    python scripts/saveExpectedBtagEff.py -i /store/cmst3/user/psilva/LJets2015/8c1e7c9/MC13TeV_TTJets${i} -o $CMSSW_BASE/src/TopLJets2015/TopAnalysis/data/expTageff${i}.root;
done
