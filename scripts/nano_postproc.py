#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
import re, time

def prepare_condor_jobs(arglist):

    inputdir = arglist[2]

    for ifile in os.listdir(inputdir):
        res = re.match("myNanoProdMc2016_NANO_([7-9]\d+).root", ifile)
        if res is None:
           continue
        fname = inputdir + "/" + ifile
        short = res.group(1)

        # prepare the log directory
        logdir = os.environ['PWD']+'/'+ "logs"
        logdir = logdir.rstrip("/")
        if not os.path.exists(logdir):
            os.system("mkdir -p "+logdir)

        # prepare the sh script
        argstring = "python"
        for iarg in range(len(arglist)):
            arg = arglist[iarg]
            if iarg == 2:
               arg += ifile
            if "--condor" in arg:
               continue
            if arglist[iarg-1] ==  "-c" or arglist[iarg-1]== "--cut":
               arg = "\"" + arg + "\""
            argstring += " " + arg
        print argstring
        condor_exec="{logdir}/{data}.sh".format(logdir=logdir, data=short)
        os.system("echo '#!/bin/bash' > {condor_exec} ; echo 'source /afs/cern.ch/work/y/yofeng/public/PFCand/CMSSW_10_2_10/src/PhysicsTools/NanoAODTools/scripts/setup_env.sh' >> {condor_exec};echo '{cmd}' >> {condor_exec} ; chmod u+x {condor_exec}".format(cmd=argstring, data=short, condor_exec=condor_exec))
       
        # write the submission script
        job_desc = """Universe = vanilla
Executable = {condor_exec}
use_x509userproxy = $ENV(X509_USER_PROXY)
Log        = {pid}.log
Output     = {pid}.out
Error      = {pid}.error
getenv      = True
environment = "LS_SUBCWD={here}"
+JobFlavour = "workday"
queue 1\n
""".format(
           here=os.environ['PWD'],
           condor_exec = condor_exec,
           pid= logdir +"/" + short,
           )
        jobdesc_file = logdir+"/"+ short+'.condor'
        print jobdesc_file
        with open(jobdesc_file,'w') as outputfile:
            outputfile.write(job_desc)
            outputfile.close()
        os.system( "condor_submit " + jobdesc_file )
        time.sleep(2)

def main():
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] outputDir inputFiles")
    parser.add_option("-s", "--postfix",dest="postfix", type="string", default=None, help="Postfix which will be appended to the file name (default: _Friend for friends, _Skim for skims)")
    parser.add_option("-J", "--json",  dest="json", type="string", default=None, help="Select events using this JSON file")
    parser.add_option("-c", "--cut",  dest="cut", type="string", default=None, help="Cut string")
    parser.add_option("-b", "--branch-selection",  dest="branchsel", type="string", default=None, help="Branch selection")
    parser.add_option("--bi", "--branch-selection-input",  dest="branchsel_in", type="string", default=None, help="Branch selection input")
    parser.add_option("--bo", "--branch-selection-output",  dest="branchsel_out", type="string", default=None, help="Branch selection output")
    parser.add_option("--friend",  dest="friend", action="store_true", default=False, help="Produce friend trees in output (current default is to produce full trees)")
    parser.add_option("--full",  dest="friend", action="store_false",  default=False, help="Produce full trees in output (this is the current default)")
    parser.add_option("--noout",  dest="noOut", action="store_true",  default=False, help="Do not produce output, just run modules")
    parser.add_option("-P", "--prefetch",  dest="prefetch", action="store_true",  default=False, help="Prefetch input files locally instead of accessing them via xrootd")
    parser.add_option("--long-term-cache",  dest="longTermCache", action="store_true",  default=False, help="Keep prefetched files across runs instead of deleting them at the end")
    parser.add_option("-N", "--max-entries", dest="maxEntries", type="long",  default=None, help="Maximum number of entries to process from any single given input tree")
    parser.add_option("--first-entry", dest="firstEntry", type="long",  default=0, help="First entry to process in the three (to be used together with --max-entries)")
    parser.add_option("--justcount",   dest="justcount", default=False, action="store_true",  help="Just report the number of selected events") 
    parser.add_option("-I", "--import", dest="imports",  type="string", default=[], action="append", nargs=2, help="Import modules (python package, comma-separated list of ");
    parser.add_option("-z", "--compression",  dest="compression", type="string", default=("LZMA:9"), help="Compression: none, or (algo):(level) ")
    parser.add_option("--condor",   dest="condor", action="store_true", default=False, help="To run the jobs on condor. One file per job");

    arglist = sys.argv
    print arglist
    (options, args) = parser.parse_args()

    if options.condor:
        print "prepare the files and submit the job"
        prepare_condor_jobs( arglist )

        return 1

    if options.friend:
        if options.cut or options.json: raise RuntimeError("Can't apply JSON or cut selection when producing friends")

    if len(args) < 2 :
	 parser.print_help()
         sys.exit(1)
    outdir = args[0]; args = args[1:]

    modules = []
    defaults_to_import =[ 
                            ('PhysicsTools.NanoAODTools.postprocessing.examples.recoilModule', 'recoilModuleConstr'),
                            ('PhysicsTools.NanoAODTools.postprocessing.examples.genPFMatchingModule', 'genPFMatchingModuleConstr'),
                        ]
    for mod, names in options.imports + defaults_to_import: 
        import_module(mod)
        obj = sys.modules[mod]
        selnames = names.split(",")
        mods = dir(obj)
        for name in selnames:
            if name in mods:
                print "Loading %s from %s " % (name, mod)
                if type(getattr(obj,name)) == list:
                    for mod in getattr(obj,name):
                        modules.append( mod())
                else:
                    modules.append(getattr(obj,name)())
    if options.noOut:
        if len(modules) == 0: 
            raise RuntimeError("Running with --noout and no modules does nothing!")
    if options.branchsel!=None:
        options.branchsel_in = options.branchsel
        options.branchsel_out = options.branchsel

    # run the PostProcessor 
    p=PostProcessor(outdir,args,
            cut = options.cut,
            branchsel = options.branchsel_in,
            modules = modules,
            compression = options.compression,
            friend = options.friend,
            postfix = options.postfix,
            jsonInput = options.json,
            noOut = options.noOut,
            justcount = options.justcount,
            prefetch = options.prefetch,
            longTermCache = options.longTermCache,
            maxEntries = options.maxEntries,
            firstEntry = options.firstEntry,
            outputbranchsel = options.branchsel_out)
    p.run()

if __name__ == "__main__":
   main()

