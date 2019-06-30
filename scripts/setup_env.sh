export SCRAM_ARCH=slc7_amd64_gcc700
cd /afs/cern.ch/work/y/yofeng/public/PFCand/CMSSW_10_2_10/src/PhysicsTools/NanoAODTools/scripts
eval $(scramv1 runtime -sh);
ulimit -c 0
