#!/usr/bin/env python
# --------------------------------------------------------
# Copyright (C) 2017 PHARMACELERA S.L.
# All rights reserved.
#
# File: calculate_ROC_AUC.py
#
# @author: arudzinska
#
# --------------------------------------------------------

# --------------------------------------------------------
# Script for generating ROC/awROC (average weighted ROC)
# analysis. Generates a curve PNG file and a CSV file with
# ROC Enrichments & AUC values for a given target. Can be
# run with option for: ROC, awROC or both at the same time
# (in the latter case creates two PNG and two CSV files).
#
# ========>  Example (ROC analysis):  <========
#
#  calculate_ROC_AUC.py -rank ranking_best_conf.txt
#    -roc actives.txt (-o my_output -log yes)
#
# -roc: TXT file with actives' names/indices (one per line)
#
# ========>  Example (awROC analysis):  <========
#
#  calculate_ROC_AUC.py -rank ranking_best_conf.txt
#    -awroc act_clustered.sdf (-o my_output -log yes)
#
# -awroc: provide original DUD LIB VS 1.0 SDF with clustered actives (in
#   case of not using DUD, the actives' file needs to be alike, i.e. contain
#   '> <id>' and '> <Cluster>' fields placed *exactly* in the same manner)
#
# ==================================================
#
# -rank is the PharmScreen ranking_best_conf.txt file
# -o argument is optional: pattern for the name of the
# output files (default: 'roc_metrics/awroc_metrics').
# -log argument is optional: write 'y' or 'yes' to 
# generate the awROC curve in semi-logarithmic scale
#
# If you want to generate both standard and average weighted ROC analysis,
# simply use -roc and -awroc options together.
#
# --------------------------------------------------------
import argparse
import ranking
import file
import matplotlib.pyplot as plt

# --------------------------------------------------------
def process_actives(actives):
    """
    Only for average weighted ROC (awROC) analysis.
    
    :param: <list> original DUD SD file with clustered actives, after processing with readlines():
    :return: molecules in clusters <dict> {molecule_id:cluster no.}, cluster structure count <dict> {cluster no.:structures}
    """
    molecules_in_clusters, cluster_structure_count = {}, {}
    for idx, line in enumerate(actives):
        if line.startswith('> <id>'):
            molecule_id = actives[idx+1].strip()
            cluster = int(actives[idx+4].strip())
            molecules_in_clusters[molecule_id] = cluster
            if cluster not in cluster_structure_count:
                cluster_structure_count[cluster] = 1
            else:
                cluster_structure_count[cluster] += 1
    return molecules_in_clusters, cluster_structure_count

# --------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Script for generating ROC/awROC analysis. Generates the (aw)ROC curve (PNG) and saves the enrichment and (aw)AUC values (CSV).')
    parser.add_argument('-rank', required=True, help='screening ranking (txt)')
    parser.add_argument('-roc', help="ROC metrics: provide TXT file with actives' names/indices (one per line)")
    parser.add_argument('-awroc', help="average weighted ROC metrics: provide SDF with clustered actives")
    parser.add_argument('-o', help='optional: naming pattern of the output files (default: roc_metrics/awroc_metrics')
    parser.add_argument('-log', help="optional: create the curve in a semi-logarithmic scale (write 'y' or 'yes')")
    args = parser.parse_args()
    
    if args.roc == None and args.awroc == None:
        print "\nERROR: type of analysis was not specified. Choose an option: -roc, -awroc or both.\n"
        exit(1)
    try:
        rank = ranking.Ranking(args.rank)
    except IOError:
        print "\nERROR: No such file or directory:", args.rank, "\n"
        exit(1)
    if args.log == 'yes' or args.log == 'y':
        log = True
    else:
        log = False
    
    try:
        if file.file_check_existance(args.roc):
            with open(args.roc) as f:
                actives = f.read().splitlines()
            print "\n|------------------------------------------------------------"
            print '| ROC ANALYSIS\n|'
            print '| Total molecules:\t', rank.total_mols
            print '| Actives retrieved:\t', len(actives)
            ROC_curve, ROCE, AUC = rank.ROC_curve(actives, log)

            # Print results to a file
            if args.o == None:
                out = "roc_metrics"
            else:
                out = str(args.o)+"_ROC"
            out_csv = open(out+'.csv', 'w')
            plt.savefig(out)
            print >> out_csv, "ROCE & AUC"
            print >> out_csv, "(C) 2017 Pharmacelera\n----------------------\n"
            print >> out_csv, "Fraction [%]\tROC Enrichment"
            for percent, enrichment in zip([0.5, 1.0, 2.0, 5.0], ROCE):
                print >> out_csv, percent, "\t", enrichment
            print >> out_csv, "\nAUC"
            print >> out_csv, AUC
            plt.gcf().clear()
            print "|\n| Done. Output: "+out+".png and "+out+".csv"
            print "|------------------------------------------------------------\n"
        else:
            print "\nERROR: File with actives (ROC analysis) not found. Stopping script.\n"
            exit(1)
    except TypeError:
        # This error occurs if the user does not specify -roc option at all
        pass
    
    try:
        if file.file_check_existance(args.awroc):
            molecules_in_clusters, cluster_structure_count = process_actives(open(args.awroc, 'r').readlines())
            print "\n|------------------------------------------------------------"
            print '| awROC ANALYSIS\n|'
            print '| Total molecules:\t', rank.total_mols
            print '| Actives retrieved:\t', len(molecules_in_clusters)
            awROC_curve, awROCE, awAUC = rank.awROC_curve(molecules_in_clusters, cluster_structure_count, log)

            # Print results to a file
            if args.o == None:
                out = "awroc_metrics"
            else:
                out = str(args.o)+"_awROC"
            out_csv = open(out+'.csv', 'w')
            plt.savefig(out)
            print >> out_csv, "awROCE & awAUC"
            print >> out_csv, "(C) 2017 Pharmacelera\n----------------------\n"
            print >> out_csv, "Fraction [%]\tawROC Enrichment"
            for percent, enrichment in zip([0.5, 1.0, 2.0, 5.0], awROCE):
                print >> out_csv, percent, "\t", enrichment
            print >> out_csv, "\nawAUC"
            print >> out_csv, awAUC
            plt.gcf().clear()
            print "|\n| Done. Output: "+out+".png and "+out+".csv"
            print "|------------------------------------------------------------\n"
        else:
            print "\nERROR: File with actives (awROC analysis) not found. Stopping script.\n"
            exit(1)
    except TypeError:
        # This error occurs if the user does not specify -awroc option at all
        pass
