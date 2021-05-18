#!/usr/bin/env python
import numpy as np
import pandas as pd
import os
import matplotlib.ticker as tick
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter, LogFormatter, LogLocator, NullFormatter, FuncFormatter, AutoMinorLocator, MultipleLocator
import atlas_mpl_style as ampl
from argparse import ArgumentParser
# from ROOT import TFile, TH1F, Double
# import ROOT
from root_numpy import hist2array
from collections import OrderedDict
from root_pandas import *
from tqdm import tqdm
from pdb import set_trace

# Set ATLAS style
ampl.use_atlas_style()
ampl.set_color_cycle(pal='ATLAS')
plt.rcParams["font.size"] = 18
plt.rcParams["axes.labelsize"] = 20
plt.rcParams["xtick.major.size"] = 7
plt.rcParams["ytick.major.size"] = 7
plt.rcParams["xtick.minor.size"] = 4
plt.rcParams["ytick.minor.size"] = 4
# plt.rcParams['xtick.minor.width'] = 0
# plt.rcParams['xtick.major.pad'] = 8

# plt.rcParams["ytick.minor.visible"] = True

def getArgs():
    """Get arguments from command line."""
    parser = ArgumentParser()
    parser.add_argument('-r', '--region', action = 'store', choices = ['all_jet', 'two_jet', 'one_jet', 'zero_jet', 'VBF'], default = 'zero_jet', help = 'Region to process')
    # parser.add_argument('--skip', action = 'store_true', default = False, help = 'skip the hadd part')
    parser.add_argument('-va', '--variable', action = 'store', choices = ['bdt', 'NN'], default = 'NN', help = 'MVA variable to use')
    parser.add_argument('-t', '--transform', action = 'store_true', default = True, help = 'use the transform scores for categroization')
    parser.add_argument("--bkgedgecolor", default="black")
    parser.add_argument("--verbose", action="store_true")
    return  parser.parse_args()


def poisson_interval(k, alpha=0.05):
    """
    uses chisquared info to get the poisson interval. Uses scipy.stats 
    (imports in function). 
    """
    from scipy.stats import chi2
    a = alpha
    low, high = (chi2.ppf(a / 2, 2 * k) / 2, chi2.ppf(1 - a / 2, 2 * k + 2) / 2)
    if k == 0:
        low = 0.0
    return low, high


def getPoissonError(n, confLevel=0.682689):
    alpha = 1. - confLevel
    lower, upper = poisson_interval(n, alpha)
    return n - lower, upper - n


def log_formatting(value):
    """
    This is for use with the FuncFormatter.

    It writes 1 and 10 as plane numbers, everything else as exponential.
    """
    roundval = round(value)
    if roundval == 1:
        base = 1
        exp = ''
    elif roundval == 10:
        base = 10
        exp = ''
    else:
        base = 10
        exp = int(round(np.log10(value)))
    return r'{}$^{{\mathdefault{{ {} }} }}$'.format(base, exp)


def ratio_formatting(value, maxDigits=2):
    """
    This is for use with the FuncFormatter.

    It writes 1 and 10 as plane numbers, everything else as exponential.
    """
    value = round(value, maxDigits)
    result = str(value)
    while result.endswith("0"):
        result = result[:-1]
    if result.endswith("."):
        result = result[:-1]
    return result


def getYields(args, region, sample, sampleConfig, edges, variable):

    file_list = [f"outputs/{region}/{process}.root" for process in sampleConfig.get("process")]
    weight = sampleConfig.get('eventWeight', 'weight')
    normFactor = sampleConfig.get('normFactor', 1)

    data = pd.DataFrame()

    for d in tqdm(read_root(file_list, key='test', columns=[variable, weight, 'm_mumu'], chunksize=500000), desc=f'Loading {sample}', bar_format='{desc}: {percentage:3.0f}%|{bar:20}{r_bar}'):

        for selection in sampleConfig.get("selection"):
            d.query(selection, inplace=True)
        data = data.append(d, ignore_index=True)

    yields = np.histogram(data[variable], bins=edges, weights=data[weight])[0]*normFactor
    err = np.sqrt(np.histogram(data[variable], bins=edges, weights=data[weight]**2)[0])*normFactor

    return yields, err


#=============================================================================#
def plot_score(args, savename, region, samples, edges, region_tag, variable, variable_name, doLogY=False, ymin=None, y_scale_sf=1.4, ratio_scale=[0.7, 1.3]):

    if ymin:
        pass
    elif doLogY:
        ymin = 0.1
    elif not doLogY:
        ymin = 0.

    fig = plt.figure(figsize=(6.1, 6), dpi=100)
    ax = plt.gca()
    # gs = gridspec.GridSpec(6, 1, hspace=0.0, wspace=0.0)
    # ax = fig.add_subplot(gs[0:6])
    # ax2 = fig.add_subplot(gs[5:6], sharex=ax)
    # plt.setp(ax.xaxis.get_label(), visible=True)
    ax.tick_params(labelsize=20)

    totalBkgYields = np.array([0.] * (len(edges) - 1))
    totalBkgErr = np.array([ymin] * (len(edges) - 1))

    for sample in samples:

        sampleConfig = samples[sample]

        yields, err = getYields(args, region, sample, sampleConfig, edges, variable)
        if args.verbose:
            print("----------------------------------------")
            print(f"Sample: {sample}")
            print(yields)
            print(err)
        # totalBkgYields, totalBkgErr = draw_sample(args, ax, yields, err, edges, sample, samples[sample], totalBkgYields, totalBkgErr, ymin)

        plotType = sampleConfig.get("plotType", "bkg")

        if plotType == "bkg":

            bottom_hist = np.maximum(totalBkgYields, ymin)
            totalBkgYields += yields
            totalBkgErr = np.sqrt(totalBkgErr**2 + err**2)

            ax.fill_between(edges,
                            np.insert(totalBkgYields, len(totalBkgYields), totalBkgYields[-1]),
                            np.insert(bottom_hist, len(bottom_hist), bottom_hist[-1]),
                            step='post', zorder=0, label=sample, facecolor=sampleConfig.get("color"), edgecolor=args.bkgedgecolor)

        if plotType == "signal":

            ax.step(edges, np.insert(yields, len(yields), yields[-1]), where='post', zorder=10, 
                    linewidth=sampleConfig.get("linewidth", 2), ls=sampleConfig.get("ls", "--"), color=sampleConfig.get("color"))

            ax.plot([], [], color=sampleConfig.get("color"), ls=sampleConfig.get("ls", "--"), label=sample + (f'($\\times{sampleConfig.get("normFactor")}$)' if sampleConfig.get("normFactor") != 1 else ''), zorder=-10)

        if plotType == "data":

            ax.plot([], [], color='black', marker='o', ls='-', label=sample, zorder=20)
            ax.errorbar(0.5 * (edges[:-1] + edges[1:]),
                         yields,
                         linestyle='None',
                         marker="o",
                         color='black',
                         ecolor='black',
                         # elinewidth=3,
                         # capsize=0,
                         yerr=err,
                         xerr=0,
                         zorder=20)

    totalBkgErr_up = totalBkgYields + totalBkgErr
    totalBkgErr_do = totalBkgYields - totalBkgErr

    ax.fill_between(edges,
                     np.insert(totalBkgErr_up, len(totalBkgErr_up), totalBkgErr_up[-1]),
                     np.insert(totalBkgErr_do, len(totalBkgErr_do), totalBkgErr_do[-1]),
                     step='post',
                     zorder=1,
                     # label='Background Uncertainty',
                     facecolor="None",
                     edgecolor='black',
                     linewidth=0.0,
                     hatch='\\\\\\\\\\')
            
    ampl.set_ylabel("Events", ax=ax)
    if doLogY:
        ax.set_yscale("log")
        ax.yaxis.set_major_locator(LogLocator(base=10, numticks=15))
        ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(1.0,10.0)*0.1, numticks=15))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ymax = np.max(totalBkgErr_up)**(1.5) / ymin**(y_scale_sf - 1)
    else:
        ymax = np.max(totalBkgErr_up) * y_scale_sf + (1 - y_scale_sf) * ymin 
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(edges[0], edges[-1])
    handles, labels = ax.get_legend_handles_labels()
    print(labels)
    handles.reverse()
    labels.reverse()
    ax.legend(handles, labels, frameon=False, prop={'size': 14}, loc='upper left', bbox_to_anchor=(0.51, 0.997)) # , prop={'size': 18} , loc=(0.51, 0.660)
    ampl.draw_atlas_label(x=0.03, y=0.96, ax=ax, status='int', energy='13 TeV', lumi=139.0, desc=r"$H \to \mu\mu$")
    ax.text(0.03, 0.830, region_tag, ha="left", va="top", multialignment="left", transform=ax.transAxes, size=14)

    ampl.set_xlabel(variable_name, ax=ax)

    # print(ax.get_xticks())

    if doLogY:
        yticklabels = [log_formatting(i) for i in ax.get_yticks()]
        ax.set_yticklabels(yticklabels)
    xticklabels = [ratio_formatting(i) for i in ax.get_xticks()]
    ax.set_xticklabels(xticklabels)

    plt.subplots_adjust(top=0.94, bottom=0.15, left=0.16, right=0.94)
    # plt.tight_layout()

    if not os.path.isdir(f'plots/{args.variable}_score'):
        os.makedirs(f'plots/{args.variable}_score')
    plt.savefig(f"plots/{args.variable}_score/{savename}.pdf", format="pdf")


def main():

    args = getArgs()

    region = "two_jet" if args.region == "VBF" else args.region
    region_tag = {"zero_jet": "0", "one_jet": "1", "two_jet": "2", "VBF": "2"}

    samples = OrderedDict([
            ('VBF', {'process': ['VBF'], 'selection': ["m_mumu >= 120", "m_mumu <= 130"], 'eventWeight': 'weight', 'normFactor': 100, 'color': "orange", 'ls': '-.', 'plotType': 'signal'}),
            ('ggF', {'process': ['ggF'], 'selection': ["m_mumu >= 120", "m_mumu <= 130"], 'eventWeight': 'weight', 'normFactor': 100, 'color': "red", 'ls': '--', 'plotType': 'signal'}),
            ('Top', {'process': ['ttbar', 'stop'], 'selection': ["m_mumu >= 120", "m_mumu <= 130"], 'eventWeight': 'weight', 'normFactor': 1, 'color': "on:cyan", 'plotType': 'bkg'}),
            ('Diboson', {'process': ['diboson'], 'selection': ["m_mumu >= 120", "m_mumu <= 130"], 'eventWeight': 'weight', 'normFactor': 1, 'color': "on:blue", 'plotType': 'bkg'}),
            ('Drell-Yan', {'process': ['Z'], 'selection': ["m_mumu >= 120", "m_mumu <= 130"], 'eventWeight': 'weight', 'normFactor': 1, 'color': "on:pink", 'plotType': 'bkg'}),
            ('Data sideband', {'process': ['data_sid'], 'selection': ["(m_mumu <= 120) | (m_mumu >= 130)"], 'eventWeight': 'weight', 'normFactor': 0.273, 'color': "on:pink", 'plotType': 'data'}),
        ])

    edges = np.linspace(0, 1, 21)

    variable = f"{args.variable}_score{'_VBF' if args.region == 'VBF' else ''}{'_t' if args.transform else ''}"
    variable_name = f"${args.variable.upper()}_{{VBF}}$" if args.region == "VBF" else f"${args.variable.upper()}_{{ggF}}^{{({region_tag[region]})}}$"    

    plot_score(args, f"{args.variable.upper()}_{args.region}", region, samples, edges, f"{region_tag[region]}-jet", variable, variable_name, doLogY=True)

    return

if __name__ == '__main__':
    main()
