import ROOT
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

def hist_integral(hist,i,j):
    err = ROOT.Double()
    if i>j:
        n=0
        err=0
    else: n = hist.IntegralAndError(i,j,err)
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

class categorizer(object):

    def __init__(self, h_sig, h_bkg):

        self.h_sig = h_sig
        self.h_bkg = h_bkg

    def fit(self, bl, br, nbin, minN=5, floatB=False, earlystop=-1, pbar=False):

        if nbin == 1:

            if floatB: return [], 0

            nsig, dsig = hist_integral(self.h_sig, bl, br)
            nbkg, dbkg = hist_integral(self.h_bkg, bl, br)

            if nbkg < minN: return -1, -1

            z, u = calc_sig(nsig, nbkg, dsig, dbkg)

            return [], z

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

