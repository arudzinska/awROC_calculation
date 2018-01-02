#!/usr/bin/env python
# ----------------------------------------------------------
# Copyright (C) 2017 PHARMACELERA S.L.
# All rights reserved.
# 
# File: ranking.py
#
# @author: arudzinska
#
# ----------------------------------------------------------
# Class to process ranking_best_conf.txt and perform ROC or
# average weighted ROC analysis on it (produces a curve,
# Enrichment values at 0.5%, 1%, 2% and 5% and AUC). Contains
# also the total number of molecules in the ranking.
#
# ----------------------------------------------------------
import os
import matplotlib.pyplot as plt
from numpy import trapz

class Ranking(object):
    """
    object - path to the TXT file with ranking_best_conf
    """
    def __init__(self, ranking_file):
        self.indices = []
        self.names = []
        ext = os.path.splitext(ranking_file)[-1].lower()
        if ext == ".txt":
            for line in open(ranking_file, 'r').readlines()[3:]:
                split = line.split()
                if line == "\n":
                    pass
                else:
                    self.indices.append(split[1])
                    self.names.append(split[9])
            if self.names[0] == 'Undefined name':
                self.names = self.indices
                print "\nWARNING: No molecules' names provided in the ranking. Treating idx of the molecule as a name.\n"
            self.total_mols = len(self.indices)
        else:
            print "\nERROR: Wrong ranking format. Only TXT.\n"
            exit(1)
        self.ROCE = []
        self.awROCE = []
        self.curve = None
        self.aw_curve = None
        self.AUC = None
        self.awAUC = None

    def plot_curve(self, x_axis, y_axis, log_scale):
        """
        Plots a curve with provided X and Y axes. Optionally can plot in 
        semi-logarithmic scale. Plot's font etc. can be adjusted in this method.
        
        :param: x_axis <list of floats>, y_axis <list of floats>, log_scale <bool>
        :return: plot <matplotlib.pyplot>
        """
        plt.plot(x_axis, y_axis)
        axes = plt.gca()
        plt.gcf().subplots_adjust(bottom=0.15)
        axes.set_xlim([0,1])
        if log_scale is True:
            axes.set_yscale('log')
        else:
            axes.set_ylim([0,1])
        # to change the font size to the default, remove 'size':'20' below:
        fsfont = {'fontname':'Serif', 'size':'20'}
        plt.xlabel('False positive fraction',**fsfont)
        plt.ylabel('True positive fraction',**fsfont)
        plt.xticks(**fsfont)
        plt.yticks(**fsfont)
        return plt
    
    def _AUC(self, x_axis, y_axis):
        """
        Calculates the awAUC (area under the curve) based on function f(x) = y/x, where y and x come from the ROC calculations.
        
        :param: x_axis <list of floats>, y_axis <list of floats>
        :return: awAUC <float>
        """
        AUC = trapz(y_axis, x_axis)
        return AUC
        
    def ROC_curve(self, actives_list, log_scale):
        """
        Generates standard ROC analysis for the ranking: the ROC curve, ROCE
        (enrichments) and AUC. Optionally can plot curves in semi-logarithmic
        scale.
        
        :param: actives_list <list of strings>, log_scale <bool>
        :return: self.curve - plot <matplotlib.pyplot>, self.ROCE <list of 4 
        floats>, self.AUC <float> 
        """
        actives_total = len(actives_list)
        decoys_total = self.total_mols - actives_total
        x_axis, y_axis = [], []
        true_positives, false_positives = 0, 0
        for molec in self.names:
            if molec in actives_list:
                true_positives += 1
            else:
                false_positives += 1
            x_axis.append(float(false_positives) / decoys_total)
            y_axis.append(float(true_positives) / actives_total)
        
        if true_positives != actives_total:
            print "\nWARNING: The actives provided as an input do not correspond to the ranking molecules. Check whether the actives' and the ranking molecules' naming is the same.\n"
        
        # ----------------------------------------------------------
        #  Plotting the ROC curve
        # ----------------------------------------------------------
        self.curve = self.plot_curve(x_axis, y_axis, log_scale)
        
        # ----------------------------------------------------------
        #  Calculating the ROC Enrichments at 4 recommended fractions.
        #  x_axis already contains the number of FP divided by total
        #  number of decoys (which gives the fraction retrieved).
        # ----------------------------------------------------------
        for percent in ([0.005, 0.01, 0.02, 0.05]):
            for fract in x_axis:
                if fract >= percent:
                    break
                else:
                    continue
            y = y_axis[x_axis.index(fract)]
            self.ROCE.append(y/fract)
            
        # ----------------------------------------------------------
        #  Calculating AUC
        # ----------------------------------------------------------
        self.AUC = self._AUC(x_axis, y_axis)

        return self.curve, self.ROCE, self.AUC
    
    def awROC_curve(self, mols_in_clusters, cluster_struct_count, log_scale):
        """
        Generates average weighted ROC analysis for the ranking: the ROC curve,
        awROCE (enrichments) and awAUC. Optionally can plot curves in 
        semi-logarithmic scale. Applies only to clustered active compounds, like
        in DUD Lib 1.0.
        
        :param mols_in_clusters: <dict> key - molecule name, value - number of
        the cluster to which it belongs
        :param cluster_struct_count: <dict> key - cluster number, value - no. of
        structures in this cluster
        :param log_scale: <bool>
        :return: self.aw_curve - plot <matplotlib.pyplot>, self.awROCE <list of 4
        floats>, self.awAUC <float>
        """
        decoys_total = self.total_mols - len(mols_in_clusters)
        x_axis, y_axis = [], []
        weighted, false_positives = 0, 0
        for molec in self.names:
            if molec in mols_in_clusters:
                cluster_number = mols_in_clusters[molec]
                structures_in_cluster = cluster_struct_count[cluster_number]
                weighted += 1.0/structures_in_cluster
            else:
                false_positives += 1
            x_axis.append(float(false_positives) / decoys_total)
            y_axis.append(float(weighted) / len(cluster_struct_count))
        
        if int(round(weighted)) != len(cluster_struct_count):
            print "\nWARNING: The actives provided as an input do not correspond to the ranking molecules. Check whether the actives' and the ranking molecules' naming is the same.\n"
        
        # ----------------------------------------------------------
        #  Plotting the awROC curve
        # ----------------------------------------------------------
        self.aw_curve = self.plot_curve(x_axis, y_axis, log_scale)
        
        # ----------------------------------------------------------
        #  Calculating the awROC Enrichments at 4 recommended fractions.
        #  x_axis already contains the number of FP divided by total
        #  number of decoys (which gives the fraction retrieved).
        # ----------------------------------------------------------
        for percent in ([0.005, 0.01, 0.02, 0.05]):
            for fract in x_axis:
                if fract >= percent:
                    break
                else:
                    continue
            y = y_axis[x_axis.index(fract)]
            self.awROCE.append(y/fract)

        # ----------------------------------------------------------
        #  Calculating awAUC
        # ----------------------------------------------------------
        self.awAUC = self._AUC(x_axis, y_axis)
        
        return self.aw_curve, self.awROCE, self.awAUC
