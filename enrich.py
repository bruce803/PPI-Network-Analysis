##Usage: python tool4fdr.py <test_folder_name> <enrich_folder_name>
##example:  python General_enrich_analysis_FDR.py WGCNA_module function_list


import sys
import os
from scipy import stats

con_list = os.listdir(sys.argv[1])
go_list = os.listdir(sys.argv[2])
outfile1 = open("./" + sys.argv[1] + sys.argv[2] + "/" + "total_result.txt", "w")
N = float(23000)


def get_rid_of_redun(l1):
    l2 = list(set(l1))
    return l2


def intersect(a, b):
    a = get_rid_of_redun([str.strip() for str in a])
    b = get_rid_of_redun([str.strip() for str in b])
    return list(set(a) & set(b))


def correct_pvalues_for_multiple_testing(pvalues, correction_type="Benjamini-Hochberg"):
    """                                                                                                   
    consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1]) 
    """
    from numpy import array, empty
    pvalues = array(pvalues)
    n = float(pvalues.shape[0])
    new_pvalues = empty(n)
    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n - rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n / rank) * pvalue)
        for i in range(0, int(n) - 1):
            if new_values[i] < new_values[i + 1]:
                new_values[i + 1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
    return new_pvalues


for con in con_list:
    ll = []
    infile1 = open("./" + sys.argv[1] + "/" + con)  # change
    one = infile1.readlines()
    m = len(one)
    for lin1 in one:
        li1 = lin1.strip().split(",")
        ll.append("%s\n" % li1[0].upper())

    outfile = open("./" + sys.argv[1] + sys.argv[2] + "/" + con[:-4] + ".txt", "w")
    outfile.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ("List1", "List2", "P_value", "FDR", "Overlap_No.", "Overlap"))
    pvalue_list = []

    ##    print "ll" CHANGE IN OUTFILE PATH ABOVE
    ##    print ll[1:5]
    for go in go_list:
        infile2 = open("./" + sys.argv[2] + "/" + go)  # change
        two = infile2.readlines()
        two = map(lambda x: x.upper(), two)
        ##        print "two"
        ##        print two[1:5]
        n = len(two)
        inter_comp = intersect(ll, two)
        k = len(inter_comp)

        if k > 1:
            k = k - 1
            a = stats.hypergeom.sf(k, N, m, n)
            pvalue_list.append(a)
        else:
            pvalue_list.append(1)
    fdr = correct_pvalues_for_multiple_testing(pvalue_list, correction_type="Benjamini-Hochberg")
    lin_n = 0

    for go in go_list:
        infile2 = open("./" + sys.argv[2] + "/" + go)
        two = infile2.readlines()
        two = map(lambda x: x.upper(), two)
        ##        print "two"
        ##        print two[1:5]
        n = len(two)
        inter_comp = intersect(ll, two)
        k = len(inter_comp)

        if k > 1:
            k = k - 1
            a = stats.hypergeom.sf(k, N, m, n)
            if fdr[lin_n] < 0.1:
                ##            print "b"
                outfile.write("%s_%s\t%s_%s\t%s\t%s\t%s\t%s\n" % (con, m, go, n, a, fdr[lin_n], k + 1, inter_comp))
                outfile1.write("%s_%s\t%s_%s\t%s\t%s\t%s\t%s\n" % (con, m, go, n, a, fdr[lin_n], k + 1, inter_comp))
        lin_n += 1

##        if stats.hypergeom.pmf(k+1,N,m,n)>stats.hypergeom.pmf(k,N,m,n):
##            if float(1-stats.hypergeom.cdf(k,N,m,n))<0.05:
##                print "a"
##                outfile.write("%s\t%s\t%s\t%s\n"%(con,go,(1-stats.hypergeom.cdf(k,N,m,n)),inter_comp))
##        else:
##             a=stats.hypergeom.sf(k,N,m,n)
##             if float(a)<0.05:
##                 print "b"
##                 outfile.write("%s\t%s\t%s\t%s\n"%(con,go,1-a,inter_comp))
outfile.close()
outfile1.close()

import argparse

parser = argparse.ArgumentParser(description='merge T&F')
parser.add_argument("--merge", "-m", action="store_true", default=False)
parser.add_argument('--fdr', '-f', type=float)

if
