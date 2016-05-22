#for i in `cat LJets2015/skim/error.txt`
#do if [[ $i == *"TTJets_scaledown"* ]]
#then name=`echo $i | cut -d '_' -f 3`
#  num=`echo $i | cut -d '_' -f 4`
#  full="MC13TeV_TTJets_${name}_${num}"
  #python scripts/runLocalAnalysis.py -i /store/cmst3/user/psilva/LJets2015/8c1e7c9/MC13TeV_TTJets_${name}/MergedMiniEvents_${num} -o LJets2015/skim/${full} -q 8nh -n 8 --method TOPWidth::RunTopWidth --tag MC13TeV_TTJets_${name}
  #root -l /store/cmst3/user/psilva/LJets2015/8c1e7c9/MC13TeV_TTJets_${name}/${full}
#fi
#done
for i in `cat list.txt`; do xrdcp ${i}*.root root://eoscms//eos/cms/store/user/byates/LJets2015/skim/${i}; done
