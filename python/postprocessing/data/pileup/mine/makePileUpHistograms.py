import ROOT

tree = ROOT.TChain("Events")
tree.Add("/eos/cms/store/user/yofeng/WRecoilNanoAOD/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/NanoAOD_0705/191014_130515/0000/myNanoProdMc2016_NANO_1.root")
tree.Add("/eos/cms/store/user/yofeng/WRecoilNanoAOD/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/NanoAOD_0705/191014_130515/0000/myNanoProdMc2016_NANO_2.root")
tree.Add("/eos/cms/store/user/yofeng/WRecoilNanoAOD/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/NanoAOD_0705/191014_130515/0000/myNanoProdMc2016_NANO_3.root")
tree.Add("/eos/cms/store/user/yofeng/WRecoilNanoAOD/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/NanoAOD_0705/191014_130515/0000/myNanoProdMc2016_NANO_4.root")

pile_trueIn = ROOT.TH1D("pile_truthIn", "pile_trueIn", 100, 0, 100.)
tree.Draw("Pileup_nTrueInt >> {H}".format(H=pile_trueIn.GetName()))

ofile = ROOT.TFile("pileupWeight_DYJets.root", "RECREATE")
pile_trueIn.SetDirectory(ofile)
ofile.Write()
ofile.Close()
