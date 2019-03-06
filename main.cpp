#include <iostream>
#include <TMinuit.h>

double z[5], x[5], y[5], error_z[5];

double func(double x, double y, double* par) {
    return ((par[0] * par[0]) / (x * x) - 1) / (par[1] + par[2] * y - par[3] * y * y);
}

void fcn(Int_t &npar, Double_t* gin, Double_t& f, Double_t* par, Int_t iflag) {
    const Int_t nBins = 5;

    Double_t chiSq = 0;
    Double_t delta;
    for (Int_t i = 0; i < nBins; i++) {
        delta = (z[i] - func(x[i], y[i], par)) / error_z[i];
        chiSq += delta * delta;
    }

    f = chiSq;
}

int main() {

    // Data
    z[0] = 1;
    z[1] = 0.96;
    z[2] = 0.89;
    z[3] = 0.85;
    z[4] = 0.78;

    double error = 0.01;
    error_z[0] = error;
    error_z[1] = error;
    error_z[2] = error;
    error_z[3] = error;
    error_z[4] = error;

    x[0] = 1.5751;
    x[1] = 1.5825;
    x[2] = 1.6069;
    x[3] = 1.6339;
    x[4] = 1.6706;

    y[0] = 1.0642;
    y[1] = 0.97685;
    y[2] = 1.13168;
    y[3] = 1.128654;
    y[4] = 1.44016;

    TMinuit* gMinuit = new TMinuit(5);
    gMinuit->SetFCN(fcn);

    double argList[10];
    int    iErFlg = 0;

    argList[0] = 1;
    gMinuit->mnexcm("SET ERR", argList, 1, iErFlg);

    // Starting values and step size
    static double vstart[4] = {3, 1, 0.1, 0.01};
    static double step[4]   = {0.1, 0.1, 0.01, 0.001};
    gMinuit->mnparm(0, "a_1", vstart[0], step[0], 0, 0, iErFlg);
    gMinuit->mnparm(1, "a_2", vstart[1], step[1], 0, 0, iErFlg);
    gMinuit->mnparm(2, "a_3", vstart[2], step[2], 0, 0, iErFlg);
    gMinuit->mnparm(3, "a_4", vstart[3], step[3], 0, 0, iErFlg);

    // Run optimization
    argList[0] = 500;
    argList[1] = 1.;
    gMinuit->mnexcm("MIGRAD", argList, 2, iErFlg);

    // Display results
    double amin, edm, errDef;
    int nvpar, nparx, icstat;
    gMinuit->mnstat(amin, edm, errDef, nvpar, nparx, icstat);
    gMinuit->mnprin(3, amin);

    return 0;
}