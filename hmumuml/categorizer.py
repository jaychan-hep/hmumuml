from ROOT import *
from math import sqrt, log, ceil
from tqdm import tqdm

def calc_sig(sig, bkg,s_err,b_err):

  ntot = sig + bkg

  if(sig <= 0): return 0, 0
  if(bkg <= 0): return 0, 0

  #Counting experiment
  significance = sqrt(2*((ntot*log(ntot/bkg)) - sig))
  #significance = sig / sqrt(bkg)

  #error on significance
  numer = sqrt((log(ntot/bkg)*s_err)**2 + ((log(1+(sig/bkg)) - (sig/bkg))*b_err)**2)
  uncert = (numer/significance)

  return significance, uncert

def hist_integral(hist,i,j,k=-999,l=-999):
    err = Double()
    if i>j or k>l:
        n=0
        err=0
    elif k == -999 and l == -999: n = hist.IntegralAndError(i,j,err)
    else: n = hist.IntegralAndError(i,j,k,l,err)
    return n, err

def sum_hist_integral(integrals,xs=1):
    err, n = 0, 0
    for integral in integrals:
        n += integral[0]
        err += integral[1]**2
    return n*xs, sqrt(err)*xs

def sum_z(zs):
    sumu=0
    sumz=0
    for i in range(len(zs)):
        sumz+=zs[i][0]**2
        sumu+=(zs[i][0]*zs[i][1])**2
    sumz=sqrt(sumz)
    sumu=sqrt(sumu)/sumz if sumz != 0 else 0
    return sumz, sumu

def copy_hist(hname, h0, bl, br):

    nbin = br - bl + 1
    l_edge = h0.GetBinLowEdge(bl)
    r_edge = l_edge + h0.GetBinWidth(bl) * nbin
    h1 = TH1F(hname, hname, nbin, l_edge, r_edge)

    for i, j in enumerate(range(bl, br+1)):

        y = h0.GetBinContent(j)
        err = h0.GetBinError(j)

        h1.SetBinContent(i+1, y)
        h1.SetBinError(i+1, err)
        
    return h1

class categorizer(object):

    def __init__(self, h_sig, h_bkg):

        self.h_sig = h_sig
        self.h_bkg = h_bkg

    def smooth(self, bl, br, SorB='B', function='Epoly2'):

        if SorB == 'S': hist = self.h_sig
        elif SorB == 'B': hist = self.h_bkg

        nbin = hist.GetSize() - 2
        l_edge = hist.GetBinLowEdge(1)
        r_edge = l_edge + hist.GetBinWidth(1) * nbin

        h_merge_list = TList()

        if bl != 1:
            h_left = copy_hist('h_left', hist, 1, bl-1)
            h_merge_list.Add(h_left)

        hist_to_smooth = copy_hist('hist_to_smooth', hist, bl, br)
        smoothed_hist = fit_BDT('smoothed_hist', hist_to_smooth, function)
        h_merge_list.Add(smoothed_hist)

        if br != nbin:
            h_right =copy_hist('h_right', hist, br+1, nbin)
            h_merge_list.Add(h_right)

        if SorB == 'S':
            hname = self.h_sig.GetName()
            self.h_sig.Delete()
            self.h_sig = TH1F(hname, hname, nbin, l_edge, r_edge)
            self.h_sig.Merge(h_merge_list)

        elif SorB == 'B':
            hname = self.h_bkg.GetName()
            self.h_bkg.Delete()
            self.h_bkg = TH1F(hname, hname, nbin, l_edge, r_edge)
            self.h_bkg.Merge(h_merge_list)

    def fit(self, bl, br, nbin, minN=5, floatB=False, earlystop=-1, pbar=False):

        if nbin == 1:

            if floatB: return [], 0

            nsig, dsig = hist_integral(self.h_sig, bl, br)
            nbkg, dbkg = hist_integral(self.h_bkg, bl, br)

            if nbkg < minN: return -1, -1

            z, u = calc_sig(nsig, nbkg, dsig, dbkg)

            return [bl], z

        elif nbin > 1:

            L = int(ceil(log(nbin, 2)))
            N2 = 2**(L - 1)
            N1 = nbin - N2

            bmax, zmax, stop = -1, -1, 0
 
            for b in (range(bl, br+1) if not pbar else tqdm(range(bl, br+1))):

                b1, z1 = self.fit(bl, b-1, N1, minN=minN, floatB=floatB, earlystop=earlystop)
                if b1 == -1: continue
                b2, z2 = self.fit(b, br, N2, minN=minN, earlystop=earlystop)
                if b2 == -1: break

                z = sqrt(z1**2 + z2**2)
                if z > zmax:
                    stop = 0
                    zmax = z
                    bmax = sorted(list(set(b1 + [b] + b2)))
                else:
                    stop += 1
                if stop == earlystop: break

            return bmax, zmax

def fit_BDT(hname, hist, function='Epoly2'):

    n_events = hist.Integral()

    nbin = hist.GetSize() - 2
    l_edge = hist.GetBinLowEdge(1)
    r_edge = l_edge + hist.GetBinWidth(1) * nbin

    bdt = RooRealVar("bdt", "bdt", l_edge, r_edge)

    a = RooRealVar("a", "a", -100, 100)
    b = RooRealVar("b", "b", -100, 100)
    c = RooRealVar("c", "c", -100, 100)
    d = RooRealVar("d", "d", -100, 100)
    e = RooRealVar("e", "e", -100, 100)
    f = RooRealVar("f", "f", -100, 100)

    pdf = {}
    pdf['power'] = RooGenericPdf("pdf","PDF for BDT distribution", "@0**@1", RooArgList(bdt, a)) #power

    # EPoly Family
    pdf['Exp'] = RooGenericPdf("pdf","PDF for BDT distribution", "exp(@1*@0)", RooArgList(bdt, a))
    pdf['Epoly2'] = RooGenericPdf("pdf","PDF for BDT distribution", "exp(@1*@0+@2*@0**2)", RooArgList(bdt, a, b))
    pdf['Epoly3'] = RooGenericPdf("pdf","PDF for BDT distribution", "exp(@1*@0+@2*@0**2+@3*@0**3)", RooArgList(bdt, a, b, c))
    pdf['Epoly4'] = RooGenericPdf("pdf","PDF for BDT distribution", "exp(@1*@0+@2*@0**2+@3*@0**3+@4*@0**4)", RooArgList(bdt, a, b, c, d))
    pdf['Epoly5'] = RooGenericPdf("pdf","PDF for BDT distribution", "exp(@1*@0+@2*@0**2+@3*@0**3+@4*@0**4+@5*@0**5)", RooArgList(bdt, a, b, c, d, e))
    pdf['Epoly6'] = RooGenericPdf("pdf","PDF for BDT distribution", "exp(@1*@0+@2*@0**2+@3*@0**3+@4*@0**4+@5*@0**5+@6*@0**6)", RooArgList(bdt, a, b, c, d, e, f))

    data = RooDataHist("dh", "dh", RooArgList(bdt), hist)
    #data.Print("all")

    pdf[function].fitTo(data, RooFit.Verbose(False), RooFit.PrintLevel(-1))

    #lam.Print()
    a.Print()
    b.Print()
    c.Print()
    d.Print()
    e.Print()
    f.Print()

    frame = bdt.frame()
    data.plotOn(frame)
    pdf[function].plotOn(frame)

    #frame.Draw()

    dof = {'power': 2, 'Exp': 2, 'Epoly2': 3, 'Epoly3': 4, 'Epoly4': 5, 'Epoly5': 6, 'Epoly6': 7}
    reduced_chi_square = frame.chiSquare(dof[function])
    probability = TMath.Prob(frame.chiSquare(dof[function]) * (nbin - dof[function]), nbin - dof[function])
    print 'chi square:', reduced_chi_square
    print 'probability: ', probability

    #raw_input('Press enter to continue')

    # fill the fitted pdf into a histogram
    hfit = TH1F(hname, hname, nbin, l_edge, r_edge)
    pdf[function].fillHistogram(hfit, RooArgList(bdt), n_events)

    return hfit
