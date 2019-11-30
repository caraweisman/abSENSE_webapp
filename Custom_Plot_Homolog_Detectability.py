import numpy as np
import glob
import sys
import matplotlib.pyplot as plt
import random
from scipy.optimize import curve_fit
from scipy.stats import chi2
import inspect
import math
import re
import sys
import glob
import datetime
import os
import warnings
from cycler import cycler
from dill.source import getsource
from scipy import stats
import matplotlib as mpl

mpl.rcParams['figure.figsize'] = [11.0, 5.0]


###### Start define functions ######


## curve to fit
def func(x, a, b):
        return a*np.exp(-b*x)

def isfloat(s):
        try:
                float(s)
                return True
        except ValueError:
                return False


# function to, where possible, use maximum likelihood estimates of a and b parameter plus estimated covariance matrix to directly sample from the probability distribution of a and b (assume Gaussian with mean of max likelihood estimates and given covariance structure)
def parameter_CI_find(mla, mlb, covar):
        
        testavals = []
        testbvals = []

        if True not in np.isinf(covar):
        
                for i in range(0, 20): # 200?
                        a = np.random.multivariate_normal([mla, mlb], covar)[0]
                        b = np.random.multivariate_normal([mla, mlb], covar)[1]
                        testavals.append(a)
                        testbvals.append(b)

                if len(testavals) > 0:
                        return testavals, testbvals 
                else:
                        return 'failed'
        else:
                return 'failed'

# function to take each of the sampled a, b values and use them to sample directly from the distribution of scores taking into account the Gaussian noise (a function of distance, a, b) 
# this gives an empirical estimate of the prediction interval 
def PI_find(testavals, testbvals, currx, bitthresh):

        # sample from score distribution: Gaussian with mean a, b and noise determined by distance (currx), a, b
        PIsamples = []
        for i in range(0, len(testavals)):
                detval = func(currx, testavals[i], testbvals[i])
                estnoise = np.sqrt(testavals[i]*(1-math.exp(-1*testbvals[i]*currx))*(math.exp(-1*testbvals[i]*currx)))
                if estnoise > 0:
                        parpairvals = []
                        for j in range(0, 20): # 200?
                                PIsamples.append(detval + np.random.normal(0, estnoise))
                                parpairvals.append(detval + np.random.normal(0, estnoise))
                else:
                        PIsamples.append(detval)
        # compute mean of sample
        mean = np.mean(PIsamples)
        # compute std dev of sample
        std = np.std(PIsamples)

        # empirically determine, from sampled scores, how many are below detectability threshold 
        undetcount = 0
        for i in range(0, len(PIsamples)):
                if PIsamples[i] < bitthresh: 
                        undetcount = undetcount + 1


        # calculate this analytically from std estimate
        pval = stats.norm.cdf(bitthresh, mean, std)
        
        # calculate 99% CI 
        (lowint, highint) = stats.norm.interval(0.99, mean, std)

        return lowint, highint, pval

###### End define functions ######

def custom_make_pred_plot(gene, seqlen, ethresh, speciesorder, rawdistances, dbsize, predspeclocs, bitscores):

        # Ignore warning that sometimes happen as a result of stochastic sampling but that doesn't affect overall computation
        warnings.filterwarnings("ignore", message="invalid value encountered in sqrt")

        # make new arrays to put truncated (not using omitted species) scores, distances
        genebitscores = []
        truncdistances = []

        dblengths = [dbsize]*len(speciesorder)

        speciesorderbydist = [x for _,x in sorted(zip(rawdistances,speciesorder))]
        pred_spec_locs = [x for _,x in sorted(zip(rawdistances,predspeclocs))]
        speciesdblengths = [x for _,x in sorted(zip(rawdistances,dblengths))]
        orderedscores = [x for _,x in sorted(zip(rawdistances,bitscores))]
        
        rawdistancesbydist = sorted(rawdistances)

        predbitscores = []
        preddistances = []
        allbitscores = []
        alldistances = []
        logscores = []

        notinfitdists = []
        notinfitspecs = []
        
        ambigspecs = []
        ambigdists = []

        absentspecs = []
        absentdists = []
        
        for k in range(0, len(orderedscores)):
                infit = False
                if orderedscores[k] != 0: 
                        allbitscores.append(float(orderedscores[k]))
                        alldistances.append(rawdistancesbydist[k])
                        logscores.append(math.log(float(orderedscores[k])))
                        #if k not in pred_spec_locs:
                        if k in pred_spec_locs:
                                predbitscores.append(float(orderedscores[k]))
                                infit = True
                                preddistances.append(rawdistancesbydist[k])
                if orderedscores[k] == 'N/A':
                        ambigspecs.append(speciesorderbydist[k])
                        ambigdists.append(rawdistancesbydist[k])
                if orderedscores[k] == 0:
                        absentspecs.append(speciesorderbydist[k])
                        absentdists.append(rawdistancesbydist[k])
                if infit == False:
                        notinfitdists.append(rawdistancesbydist[k])
                        notinfitspecs.append(speciesorderbydist[k])

        #print predbitscores
        #print preddistances
        #print allbitscores
        #print alldistances

        if len(predbitscores) < 3:
                sys.exit('Too few data points! There must be at least 3. Quitting. \n') 
                                        
        smoothx = np.linspace(min(rawdistances), max(rawdistances), 1000)
        predictions = []
        highs = []
        lows = []

        notinfitpvals = []

        # define average database size for general line
        #avgdbsize = np.mean(speciesdblengths.astype(float))
        avgdbsize = np.mean(speciesdblengths)
        globalbitthresh = -1*math.log(ethresh/(avgdbsize * seqlen), 2)
        bitthresh = globalbitthresh

        try: 
                (a, b), covar = curve_fit(func, preddistances, predbitscores)
                slope, intercept, r_value, p_value, std_err = stats.linregress(alldistances,logscores)
        except RuntimeError:
                pass
##                print 'Runtime Error'
        parout = parameter_CI_find(a, b, covar) 
        if parout != 'failed':
                testavals, testbvals = parout
                for j in range(0, len(smoothx)):
                        predictions.append(round(func(smoothx[j], a,b),2))
                        lowprediction, highprediction, pval = PI_find(testavals, testbvals, smoothx[j], bitthresh)
                        highs.append(highprediction)
                        lows.append(lowprediction)
##                print '\n'
##                print '\n'
##                print '\n'
##                print 'Results:'
##
##                print '\n'
##                print '\n'
##                print 'Best-fit a parameter from bitscores included in fit (black points): ', round(a,2)
##                print 'Best-fit b parameter from bitscores included in fit (black points): ', round(b,2)
##                print 'r squared from all bitscores (black points and orange points, if any):', round(r_value**2,2)
##                print '\n'
##                print '\n'
                
                aparam = str(round(a,2))
                bparam = str(round(b, 2))
                rsq = str(round(r_value**2,2))
                mlpreds = []
                highnns = []
                lownns = []
                undets = []
                for j in range(0, len(notinfitspecs)):
                        dblen = float(speciesdblengths[speciesorderbydist.index(notinfitspecs[j])])
                        bitthresh = -1*math.log(ethresh/(dblen * seqlen), 2)
                        lowprediction, highprediction, pval = PI_find(testavals, testbvals, notinfitdists[j], bitthresh)
                        #print notinfitspecs[j], ':'
                        #print 'Maximum likelihood bitscore prediction: ', str(round(func(notinfitdists[j], a,b),2))
                        mlpreds.append(str(round(func(notinfitdists[j], a,b),2)))
                        #print '99th percentile (high) prediction: ', str(round(highprediction, 2))
                        highnns.append(str(round(highprediction, 2)))
                        #print '1st percentile (low) prediction: ', str(round(lowprediction, 2))
                        lownns.append(str(round(lowprediction, 2)))
                        #print 'Probability of homolog being undetected: ', pval
                        undets.append(str(round(pval, 2)))
                        #print '\n'

        startheight = 0.88

        labels = speciesorderbydist
        afont = {'fontname':'Arial'}
        plt.title('Gene: ' + gene + '\n' + 'a = ' + str(round(a,1)) + ', b = ' + str(round(b,2)) + '\n' + '$r^2$ = ' + str(round(r_value**2, 2)), color='black', fontsize=13, fontweight='bold', **afont)
        plt.scatter(preddistances, predbitscores, s=40, c='black', label='Bitscores of detected orthologs')
        plt.plot(smoothx, predictions, c='red', label='Predicted bitscore')
        plt.plot(smoothx, highs, c='black')
        plt.plot(smoothx, lows, c='black')
        plt.fill_between(smoothx, highs, lows, facecolor='blue', alpha=0.2, label='99% confidence interval')
        plt.ylabel('Bitscore', fontsize=13, labelpad=10, **afont)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().tick_params(axis='x', width=2, length=7, direction='inout')
        plt.xticks(rawdistancesbydist, labels, fontsize=10, rotation=90) #labels
        totrange = max(smoothx) - min(smoothx)
        inc = float(totrange)/100
        for i in range(0, len(absentspecs)):
                if i == 0:
                        plt.axvspan(absentdists[i] - inc, absentdists[i] + inc, facecolor='#fc8123', alpha=0.3, label='No homolog detected in species', capstyle='round')
                        plt.gca().get_xticklabels()[labels.index(absentspecs[i])].set_color('#fc8123')
                        plt.gca().get_xticklabels()[labels.index(absentspecs[i])].set_weight('bold')
                        
                else:
                        plt.axvspan(absentdists[i] - inc, absentdists[i] + inc, facecolor='#fc8123', alpha=0.3, capstyle='round')
                        plt.gca().get_xticklabels()[labels.index(absentspecs[i])].set_color('#fc8123')
                        plt.gca().get_xticklabels()[labels.index(absentspecs[i])].set_weight('bold')
        for i in range(0, len(ambigspecs)):
                if i == 0:
                        plt.gca().get_xticklabels()[labels.index(ambigspecs[i])].set_color('#a3a29b')
                        
                else:
                        plt.gca().get_xticklabels()[labels.index(ambigspecs[i])].set_color('#a3a29b')
        plt.yticks(fontsize=10)
        #plt.figtext(0.5,startheight, 'Gene = ' + gene, color='black', fontsize=12, fontweight='bold', **afont)
        #plt.figtext(0.5,startheight - 0.05,'a = ' + str(round(a,1)) + ', b = ' + str(round(b,2)) )
        #plt.figtext(0.5,startheight - 0.1,' $r^2$ = ' + str(round(r_value**2, 2)))
        plt.axhline(y=globalbitthresh, linestyle='dashed', c='black', label='Detectability threshold')
        plt.xlim([-inc, max(rawdistances)+inc])
        plt.ylim([0, max(predictions)+max(predictions)/10])
        #print gene
        #print absentspecs
        handles, labels = plt.gca().get_legend_handles_labels()
        if len(absentspecs) > 0: 
                order = [3, 0, 4, 1, 2]
        else:
                order = [2,0, 3,1]
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], fontsize=9)

        #plt.show()
        #plt.savefig('testoutinscript.png')
        return plt, aparam, bparam, rsq, notinfitspecs, mlpreds, highnns, lownns, undets, ambigspecs, absentspecs, speciesorderbydist, orderedscores




#make_pred_plot()
#plt.show()






