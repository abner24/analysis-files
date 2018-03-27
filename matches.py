#!/bin/env python3
# ##################################################
# Title : re_pe.py
# Author : Lim Ming Sun (Abner)
# Date : Feb 05 2018 16:53:50
# Description : python script to calculate precision and recall rate
# ##################################################

import os
import io
import sys
import gzip
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib import axes
from copy import copy
#from StringIO import StringIO

#original VCF
#vcf_ori=sys.argv[1]
# the snv file from phylovar
snv_file=sys.argv[1]
# VCF file
vcf_file=sys.argv[2]
out_dir=os.path.dirname(vcf_file)

class variants:
    def __init__(self):
        self.chr = []
        self.pos = []
        self.freq = []

def main():
        
        snv = variants()
        with open(snv_file,'r') as f_read:
                for line in f_read:
                        if line[0] !='#':
                                snv.chr.append(line.split('\t')[0])
                                snv.pos.append(int(line.split('\t')[1]))
                                snv.freq.append(float(line.split('\t')[4]))
        f_read.close()
        # to change from snv from base 0 to 1.
        # snv.pos=np.array(snv.pos)+np.ones(len(pos),dtype=int)
        
        vcf = variants()
        with io.TextIOWrapper(io.BufferedReader(gzip.open(vcf_file,'rb'))) as vcf_read:
                for line in vcf_read:
                        if line[0]!='#':
                                vcf.chr.append(line.split('\t')[0])
                                vcf.pos.append(int(line.split('\t')[1]))
                                vcf.freq.append(float(line.split('\t')[10].split(':')[4]))
                                #print(line.split(b'\t')[1])
        vcf_read.close()

        match=0
        #pos_vcf=np.array(pos_vcf)
        out_file=os.path.join(out_dir,'unaccounted.txt')
        undetected_freq=copy(snv.freq)
        detected_freq=[]
        vcf_detected=[]
        vcf_undetected=copy(vcf.freq)
        snv_tracker=[]
        vcf_tracker=[]
        with open(out_file,'w') as f:
                f.write('chromosome and location of of SNV unaccounted for\n')
                f.write('#chrom\tlocation\tfreq\t\n')
                for i in range(len(snv.pos)):
                    for j in range(len(vcf.pos)):
                        if (snv.pos[i]+1) == vcf.pos[j] and snv.chr[i] == vcf.chr[j]:
                            detected_freq.append(snv.freq[i])
                            vcf_detected.append(vcf.freq[j])
                            vcf_undetected.remove(vcf.freq[j])
                            undetected_freq.remove(snv.freq[i])
                            match+=1
                            f.write(str(snv.chr[i])+'\t'+str(snv.pos[i])+'\t'+str(snv.freq[i]))
                            break
        f.close()

        recall = match/len(snv.pos)
        precision = match/len(vcf.pos)
        print('number of matches found '+str(match))
        print('the recall is ' +str(recall))
        print('the precision is '+str(precision))
        print('there are {} false positives'.format(len(vcf_undetected)))
        print('there are {} variants in phylovar not recalled'.format(len(undetected_freq)))
        return detected_freq, undetected_freq, vcf_detected, vcf_undetected;



y0,y1,y2,y3=main()

plt.figure(1)
freq,bins_freq,patches=plt.hist(y0,10)
ax1=plt.gca()
ax1.axes.tick_params('x',length=0.05)
ax1.axes.tick_params('y',length=1)
plt.xlabel('frequency')
plt.ylabel('counts')
plt.title('SNVs detected by mutect simulated by CSiTE')
plt.savefig('snvs.png')
plt.figure(2)
unde_freq,bins_unde,patches_unde=plt.hist(y1,10)
ax2=plt.gca()
ax2.axes.tick_params('x',length=0.05)
ax1.axes.tick_params('y',length=1)
plt.xlabel('freq')
plt.ylabel('counts')
plt.title('SNVs undetected by mutect simulated by CSiTE')
plt.savefig('snvs_undetected.png')
plt.figure(3)
freq,bins_freq,patches=plt.hist(y2,10)
ax3=plt.gca()
ax3.axes.tick_params('x',length=0.05)
ax3.axes.tick_params('y',length=1)
plt.xlabel('frequency')
plt.ylabel('counts')
plt.title('Mutect PASS filtered output detected')
plt.savefig('vcf_detected.png')
plt.figure(4)
freq,bins_freq,patches=plt.hist(y3,bins="auto")
ax4=plt.gca()
ax4.axes.tick_params('x',length=0.05)
ax4.axes.tick_params('y',length=1)
plt.xlabel('frequency')
plt.ylabel('counts')
plt.title('Mutect PASS filtered output undetected')
plt.savefig('vcf_undetected.png')
