#ifndef HiggsAnalysis_CombinedLimit_PdfDiagonalizer_h
#define HiggsAnalysis_CombinedLimit_PdfDiagonalizer_h

#include <RooAbsPdf.h>
#include <RooAddition.h>
#include <RooCustomizer.h>
#include <RooFitResult.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <cstdio>

// class RooWorkspace;
// struct RooFitResult;
// struct RooAbsPdf;

#include <RooArgList.h>
#include <RooCmdArg.h>
#include <string>

class PdfDiagonalizer {
public:
  PdfDiagonalizer();
  PdfDiagonalizer(const char *name, RooWorkspace *w, RooFitResult &result);

  RooAbsPdf *diagonalize(RooAbsPdf &pdf);
  RooAbsPdf *diagonalizeWithEigenVariations(RooAbsPdf &pdf,
                                            RooFitResult &result, int index_eig,
                                            int nSigma);
  const RooArgList &originalParams() { return parameters_; }
  const RooArgList &diagonalParams() { return eigenVars_; }

private:
  std::string name_;
  RooArgList parameters_;
  RooArgList eigenVars_;
  RooArgList replacements_;
  // TVectorD& eigen_gl;
};

#endif
