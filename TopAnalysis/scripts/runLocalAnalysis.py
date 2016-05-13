import os
import sys
import optparse
import ROOT
import pickle
import json
from TopLJets2015.TopAnalysis.storeTools import *

"""
Wrapper to be used when run in parallel
"""
def RunMethodPacked(args):

    method,inF,outF,channel,charge,wgtH,flav,runSysts=args
    method=method if method.find('::')<0 else method.split('::')[1]
    print 'Running ',method,' on ',inF
    print 'Output file',outF
    print 'Selection ch=',channel,' charge=',charge,' flavSplit=',flav,' systs=',runSysts
    if wgtH : print 'Weight histogram is available'

    try:
        getattr(ROOT,method)(str(inF),str(outF),channel,charge,flav,wgtH,runSysts)
    except :
        print 50*'<'
        print "  Problem  (%s) with %s continuing without"%(sys.exc_info()[1],inF)
        print 50*'<'
        return False
    return True

"""
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-m', '--method',      dest='method',      help='method to run [%default]',                   default='TOP-16-006::RunTop16006',  type='string')
    parser.add_option('-i', '--in',          dest='input',       help='input directory with files or single file [%default]',  default=None,       type='string')
    parser.add_option('-o', '--out',         dest='output',      help='output directory (or file if single file to process)  [%default]',  default='analysis', type='string')
    parser.add_option(      '--only',        dest='only',        help='csv list of samples to process  [%default]',             default=None,       type='string')
    parser.add_option(      '--runSysts',    dest='runSysts',    help='run systematics  [%default]',                            default=False,      action='store_true')
    parser.add_option(      '--cache',       dest='cache',       help='use this cache file  [%default]',                        default='data/genweights.pck', type='string')
    parser.add_option(      '--flav',        dest='flav',        help='split according to heavy flavour content  [%default]',   default=0,          type=int)
    parser.add_option(      '--ch',          dest='channel',     help='channel  [%default]',                                    default=13,         type=int)
    parser.add_option(      '--charge',      dest='charge',      help='charge  [%default]',                                     default=0,          type=int)
    parser.add_option(      '--tag',         dest='tag',         help='normalize from this tag  [%default]',                    default=None,       type='string')
    parser.add_option('-q', '--queue',       dest='queue',       help='submit to this queue  [%default]',                       default='local',    type='string')
    parser.add_option('-n', '--njobs',       dest='njobs',       help='# jobs to run in parallel  [%default]',                                default=0,    type='int')
    (opt, args) = parser.parse_args()

    #compile macro
    ROOT.FWLiteEnabler.enable()
    ROOT.gSystem.Load('libTopLJets2015TopAnalysis.so')
    srcCode=opt.method.split('::')[0]
    ROOT.gROOT.LoadMacro('src/%s.cc+'%srcCode)
    
    #parse selection list
    onlyList=[]
    try:
        onlyList=opt.only.split(',')
    except:
        pass

    #prepare output if a directory
    if not '.root' in opt.output:
        os.system('mkdir -p %s'%opt.output)

    #read normalization
    cachefile = open(opt.cache, 'r')
    genWgts   = pickle.load(cachefile)
    cachefile.close()        
    print 'Normalization read from cache (%s)' % opt.cache
    
    #process tasks
    task_list = []
    processedTags=[]
    if '.root' in opt.input:
        inF=opt.input
        if '/store/' in inF and not 'root:' in inF : inF='root://eoscms//eos/cms'+opt.input        
        outF=opt.output
        wgtH=None
        if opt.tag :
            if opt.tag in genWgts:
                wgtH=genWgts[opt.tag]
        print inF,outF,opt.channel,opt.charge,wgtH,opt.flav,opt.runSysts
        task_list.append( (opt.method,inF,outF,opt.channel,opt.charge,wgtH,opt.flav,opt.runSysts) )
    else:

        inputTags=getEOSlslist(directory=opt.input,prepend='')
        for baseDir in inputTags:

            tag=os.path.basename(baseDir)
            if tag=='backup' : continue

            #filter tags
            if len(onlyList)>0:
                processThisTag=False
                for itag in onlyList:
                    if itag in tag:
                        processThisTag=True
                if not processThisTag : continue

            wgtH=genWgts[opt.tag] if opt.queue=='local' else tag
            input_list=getEOSlslist(directory='%s/%s' % (opt.input,tag) )
            for ifile in xrange(0,len(input_list)):
                inF=input_list[ifile]
                outF=os.path.join(opt.output,'%s_%d.root' %(tag,ifile))
                #doFlavourSplitting=True if ('MC13TeV_WJets' in tag or 'MC13TeV_DY50toInf' in tag) else False
                #if doFlavourSplitting:
                #    for flav in [0,1,4,5]:
                #        task_list.append( (method,inF,outF,opt.channel,opt.charge,wgtH,flav,opt.runSysts) )
                #else:
                task_list.append( (opt.method,inF,outF,opt.channel,opt.charge,wgtH,0,opt.runSysts) )

    #run the analysis jobs
    if opt.queue=='local':
        print 'launching %d tasks in %d parallel jobs'%(len(task_list),opt.njobs)
        if opt.njobs == 0:
            for args in task_list: RunMethodPacked(args)
        else:
            from multiprocessing import Pool
            pool = Pool(opt.njobs)
            pool.map(RunMethodPacked, task_list)
    else:
        print 'launching %d tasks to submit to the %s queue'%(len(task_list),opt.queue)
        cmsswBase=os.environ['CMSSW_BASE']
        for method,inF,outF,channel,charge,tag,flav,runSysts in task_list:
            localRun='python %s/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py -i %s -o %s --charge %d --ch %d --tag %s --flav %d --method %s' % (cmsswBase,inF,outF,charge,channel,tag,flav,method)
            if runSysts : localRun += ' --runSysts'            
            cmd='bsub -q %s %s/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh \"%s\"' % (opt.queue,cmsswBase,localRun)
            os.system(cmd)

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
