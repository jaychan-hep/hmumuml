from ROOT import *
from math import sqrt, log
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-r', '--region', action = 'store', choices = ['all_jet', 'two_jet', 'one_jet', 'zero_jet'], default = 'zero_jet', help = 'Region to process')
args = parser.parse_args()

def calc_sig(sig, bkg):

  ntot = sig + bkg
  significance = sqrt(2*((ntot*log(ntot/bkg)) - sig))

  return significance


gROOT.SetBatch(True)

region = args.region

f_sig = TFile('outputs_2D/model_%s/sig.root'%region)
f_VBF = TFile('outputs_2D/model_%s/VBF.root'%region)
f_ggF = TFile('outputs_2D/model_%s/ggF.root'%region)
f_data = TFile('outputs_2D/model_%s/data.root'%region)
t_sig = f_sig.Get('test')
t_VBF = f_VBF.Get('test')
t_ggF = f_ggF.Get('test')
t_data = f_data.Get('test')

h_sig = TH1F('h_sig','h_sig',18,1,19)
h_VBF = TH1F('h_VBF','h_VBF',18,1,19)
h_ggF = TH1F('h_ggF','h_ggF',18,1,19)
h_data = TH1F('h_data','h_data',18,1,19)

t_sig.Draw('bdt_category>>h_sig','weight*(m_mumu>=120&&m_mumu<=130)')
t_VBF.Draw('bdt_category>>h_VBF','weight*(m_mumu>=120&&m_mumu<=130)')
t_ggF.Draw('bdt_category>>h_ggF','weight*(m_mumu>=120&&m_mumu<=130)')
t_data.Draw('bdt_category>>h_data','weight*(235005./769468.)*((m_mumu>=110&&m_mumu<=160)&&!(m_mumu>=120&&m_mumu<=130))')

print 'Total VBF events:    %.3f'%h_VBF.Integral()
print 'Total ggF events:    %.3f'%h_ggF.Integral()
print 'Total signal events: %.3f'%h_sig.Integral()
print 'Total BKG events: %.3f'%h_data.Integral()


print 'CAT NO.  VBF Events  ggF Events  Signal Events  BKG Events  VBF purity  S/B ratio  Significance'
for i in range(1,11 if region == 'two_jet' else 6):
    purity = h_VBF.Integral(i,i)/h_sig.Integral(i,i)*100
    sbratio = h_sig.Integral(i,i)/h_data.Integral(i,i)
    z = calc_sig(h_sig.Integral(i,i), h_data.Integral(i,i)) 
    print 'CAT%d:       %.3f         %.3f          %.3f        %.3f        %.3f%%        %.4f      %.2f'%(i, h_VBF.Integral(i,i), h_ggF.Integral(i,i), h_sig.Integral(i,i), h_data.Integral(i,i), purity, sbratio, z)
