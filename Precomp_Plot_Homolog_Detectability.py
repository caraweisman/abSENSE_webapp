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
import matplotlib.style
import matplotlib as mpl



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
def parameter_CI_find(mla, mlb, covar, samplesize):
        
        testavals = []
        testbvals = []

        if True not in np.isinf(covar):
        
                for i in range(0, samplesize): # 200?
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
def PI_find(testavals, testbvals, currx, bitthresh, samplesize):

        # sample from score distribution: Gaussian with mean a, b and noise determined by distance (currx), a, b
        PIsamples = []
        for i in range(0, len(testavals)):
                detval = func(currx, testavals[i], testbvals[i])
                estnoise = np.sqrt(testavals[i]*(1-math.exp(-1*testbvals[i]*currx))*(math.exp(-1*testbvals[i]*currx)))
                if estnoise > 0:
                        parpairvals = []
                        for j in range(0, samplesize): # 200?
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

        # compute fraction of sampled scores below threshold = P(undetected) = empriical "p-value"
        emppval = float(undetcount)/float(len(PIsamples))

        # calculate this analytically from std estimate
        pval = stats.norm.cdf(bitthresh, mean, std)
        
        # calculate 99% CI 
        (lowint, highint) = stats.norm.interval(0.99, mean, std)

        return lowint, highint, pval, emppval

###### End define functions ######

def fungi_make_pred_plot(gene, pred_specs, clade):

        if clade == 'fungi':
                distancefile = np.genfromtxt('Fungi_Data/Fungi_Distances', dtype=str, delimiter='\t')
                bitscores = np.genfromtxt('Fungi_Data/Fungi_Bitscores', dtype=str, delimiter='\t')
                genelengths = np.genfromtxt('Fungi_Data/S_cer_Protein_Lengths', dtype = str, delimiter='\t')
                speciesdblengths = np.genfromtxt('Fungi_Data/Fungi_Database_Lengths', dtype = str, delimiter='\t')
                mpl.rcParams['figure.figsize'] = [9.0, 5.0]

        elif clade == 'insects':

                distancefile = np.genfromtxt('Insect_Data/Insect_Distances', dtype=str, delimiter='\t')
                bitscores = np.genfromtxt('Insect_Data/Insect_Bitscores', dtype=str, delimiter='\t')
                genelengths = np.genfromtxt('Insect_Data/D_mel_Protein_Lengths', dtype = str, delimiter='\t')
                speciesdblengths = np.genfromtxt('Insect_Data/Insect_Database_Lengths', dtype = str, delimiter='\t')
                mpl.rcParams['figure.figsize'] = [11.0, 5.0]

        ethresh = 0.001

        ## Get distances, species from distance file; gene list from bitscore file
        speciesorder = distancefile[0]
        rawdistances = distancefile[1].astype(float)
        genelist = bitscores[1:,0] # skip header

        ## Ensure that species order and distances are given in ascending distance order (not important for computation or text output, but important for visualization)
        ## Also determine the locations in the ordering to be used of the species to omit from curve fit
                                
        speciesorderbydist = [x for _,x in sorted(zip(rawdistances,speciesorder))]
        rawdistancesbydist = sorted(rawdistances)

        speciestotallengths = []
        for i in range(0, len(speciesorderbydist)):
                for j in range(0, len(speciesdblengths)):
                        if speciesorderbydist[i] in speciesdblengths[j][0]:
                                speciestotallengths.append(float(speciesdblengths[j][1]))

        pred_spec_locs = []
        for i in range(0, len(speciesorderbydist)): 
                if speciesorder[i] in pred_specs:
                        pred_spec_locs.append(i)

        # Ignore warning that sometimes happen as a result of stochastic sampling but that doesn't affect overall computation
        warnings.filterwarnings("ignore", message="invalid value encountered in sqrt")

        # make new arrays to put truncated (not using omitted species) scores, distances
        genebitscores = []
        truncdistances = []

        # use value given above as gene length

        for i in range(0, len(genelengths)):
                if gene in genelengths[i][0]:
                        seqlen = float(genelengths[i][1])
        

        speciesorderbydist = [x for _,x in sorted(zip(rawdistances,speciesorder))]
        rawdistancesbydist = sorted(rawdistances)

        predbitscores = []
        preddistances = []
        predspecs = []
        
        allbitscores = []
        alldistances = []
        logscores = []

        notinfitdists = []
        notinfitspecs = []

        ambigspecs = []
        ambigdists = []

        absentspecs = []
        absentdists = []
        for j in range(0, len(bitscores)): 
                if gene in bitscores[j][0]:
                        orderedscores = [x for _,x in sorted(zip(rawdistances,bitscores[j][1:]))]
                        #print orderedscores
                        for k in range(0, len(orderedscores)):
                                infit = False
                                if isfloat(orderedscores[k]) == True and orderedscores[k] != '0': 
                                        allbitscores.append(float(orderedscores[k]))
                                        alldistances.append(rawdistancesbydist[k])
                                        logscores.append(math.log(float(orderedscores[k])))
                                        #if k not in pred_spec_locs:
                                        if k in pred_spec_locs:
                                                predbitscores.append(float(orderedscores[k]))
                                                infit = True
                                                preddistances.append(rawdistancesbydist[k])
                                                predspecs.append(speciesorderbydist[k])
                                if orderedscores[k] == 'N/A':
                                        ambigspecs.append(speciesorderbydist[k])
                                        ambigdists.append(rawdistancesbydist[k])
                                if orderedscores[k] == '0':
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
        avgdbsize = np.mean(speciesdblengths[:,1].astype(float))
        globalbitthresh = -1*math.log(ethresh/(avgdbsize), 2)
        bitthresh = globalbitthresh

        try: 
                (a, b), covar = curve_fit(func, preddistances, predbitscores)
                slope, intercept, r_value, p_value, std_err = stats.linregress(alldistances,logscores)
        except RuntimeError:
                pass
                #print 'Runtime Error'
        parout = parameter_CI_find(a, b, covar, 20) 
        if parout != 'failed':
                testavals, testbvals = parout
                for j in range(0, len(smoothx)):
                        predictions.append(round(func(smoothx[j], a,b),2))
                        lowprediction, highprediction, pval, emppval = PI_find(testavals, testbvals, smoothx[j], bitthresh, 20)
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
                for j in range(0, len(absentdists)):
                        print (bitthresh)
                        for k in range(1, len(speciesdblengths)):
                                # now use species-specific db length
                                if absentspecs[j] in speciesdblengths[k][0]:
                                        dblen = float(speciesdblengths[k][1])
                                        bitthresh = -1*math.log(ethresh/(dblen), 2)
                        lowprediction, highprediction, pval, emppval = PI_find(testavals, testbvals, absentdists[j], bitthresh, 20)
                        print (bitthresh, pval, emppval)
                        #print absentspecs[j], ':'
                        #print 'Maximum likelihood bitscore prediction: ', str(round(func(absentdists[j], a,b),2))
                        mlpreds.append(str(round(func(absentdists[j], a,b),2)))
                        #print '99th percentile (high) prediction: ', str(round(highprediction, 2))
                        highnns.append(str(round(highprediction, 2)))
                        #print '1st percentile (low) prediction: ', str(round(lowprediction, 2))
                        lownns.append(str(round(lowprediction, 2)))
                        #print 'Probability of homolog being undetected: ', str(round(pval, 2))
                        undets.append(str(round(pval, 2)))
                        #print '\n'


        startheight = 0.88

        labels = speciesorderbydist
##        plt.annotate(str(rawdistancesbydist[0]), (rawdistancesbydist[0], 0))
##        plt.annotate(str(rawdistancesbydist[len(rawdistancesbydist)-1]), (rawdistancesbydist[len(rawdistancesbydist)-1], 0))

##        labels = []
##        for k in range(0, len(speciesorderbydist)):
##                if k == 0 or k == len(speciesorderbydist)-1:
##                        labels.append(speciesorderbydist[k] + ' (' + str(rawdistancesbydist[k]) + ')')
##                else:
##                        labels.append(speciesorderbydist[k])
##        print labels
        afont = {'fontname':'Arial'}
        #if len(alldistances) > len(preddistances): 
                #plt.scatter(alldistances, allbitscores, s=40, c='orange', label='Bitscores available but not used in fit')
        plt.title('Gene: ' + gene + '\n' + 'a = ' + str(round(a,1)) + ', b = ' + str(round(b,2)) + '\n' + '$r^2$ = ' + str(round(r_value**2, 2)), color='black', fontsize=13, fontweight='bold', **afont)
        plt.scatter(preddistances, predbitscores, s=40, c='black', label='Bitscores of detected orthologs')
        plt.plot(smoothx, predictions, c='red', label='Predicted bitscore')
        plt.plot(smoothx, highs, c='black')
        plt.plot(smoothx, lows, c='black')
        plt.fill_between(smoothx, highs, lows, facecolor='blue', alpha=0.2, label='99% confidence interval')
        if clade == 'fungi':
                plt.xlabel('Evolutionary distance from S. cerevisiae (substitutions/site)', fontsize=13, **afont)
        elif clade == 'insects':
                plt.xlabel('Evolutionary distance from D. melanogaster (substitutions/site)', fontsize=13, **afont)
        plt.ylabel('Bitscore', fontsize=13, labelpad=10, **afont)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().tick_params(axis='x', width=2, length=7, direction='inout')
        plt.xticks(rawdistancesbydist, labels, fontsize=10, rotation=90) #labels
        for i in range(0, len(absentspecs)):
                if i == 0:
                        plt.axvspan(absentdists[i] - 0.01, absentdists[i] + 0.01, facecolor='#fc8123', alpha=0.3, label='No homolog detected in species', capstyle='round')
                        plt.gca().get_xticklabels()[labels.index(absentspecs[i])].set_color('#fc8123')
                        plt.gca().get_xticklabels()[labels.index(absentspecs[i])].set_weight('bold')
                        
                else:
                        plt.axvspan(absentdists[i] - 0.01, absentdists[i] + 0.01, facecolor='#fc8123', alpha=0.3, capstyle='round')
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
        plt.axhline(y=globalbitthresh, linestyle='dashed', c='black', label='Detectability threshold (E=0.001)')
        plt.xlim([-0.02, max(rawdistances)+0.04])
        plt.ylim([0, max(predictions)+max(predictions)/5])
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






