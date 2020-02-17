#!/net/fantasia/home/sarahgra/Python-3.8.0/python

import timeit
import argparse
from subprocess import Popen, PIPE
import gzip
#import io
import numpy as np
from glob import glob

#This python script can be used to calculate genetic risk scores from UKB dosages in VCF file
#This is only intended for chr1-22, as it does not account for any differences in dosage handling for males vs. females
#Be careful as this does not check the order of samples in the vcf files and assumes that they are all the same.  This will mostly be an issue for X chromosome vcf files compared to autosome.  
parser = argparse.ArgumentParser(description='Enter files to use for PRS calculation')
parser.add_argument('-s', '--score_file')
parser.add_argument('-t', '--test_dosages', nargs='*', help="Enter a full path to the UKB vcf files, eg. *_v3_s486743.vcf.gz")
parser.add_argument("-i", '--id_vcf', help="Enter a full path to the UKB vcf file to be used to extract sample IDs, eg. chr1.vcf.gz")
parser.add_argument('-o', '--output_file', default="polygenic_risk_score_results.txt", help="Enter the file name for results")
args=parser.parse_args()


#Create dictionary of scores per variant
#dictionary format variant:(effectallele, effect)
with open(args.score_file) as f:
    score_file_dict = {line.split()[0]:(line.split()[1],float(line.split()[2])) for line in f}

#Create dictionary to keep track of total scores per person, set initial value to zero
sample_id = Popen(['bcftools', 'query',  '-l', args.id_vcf], stdout=PIPE, universal_newlines=True)
sample_id = sample_id.communicate()[0].rstrip().split("\n")
sample_score_dict = {x:0 for x in sample_id}    

#flatten nested lists
def flatten(l, ltypes=(list, tuple)):
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)

#Write out regions file for tabix 
with open(args.score_file) as f:
    regions = np.genfromtxt(f, delimiter=":", usecols=(0,1), names=("Chr", "Pos"), dtype=None)

regions = np.sort(regions, order=["Chr", "Pos"])

regions_output_name = "Regions_" + args.score_file
with open(regions_output_name, 'w') as out:
    np.savetxt(regions_output_name, regions, delimiter="\t", fmt='%d')
        
#Calculate dosages per individual
#Make sure to check allele
file_list = glob(args.test_dosages)

start_time = timeit.default_timer()

for file_x in file_list:
    print("Now reading in: %s" % file_x)
    #bcftools will pull out the dosage directly
    #f = Popen(['tabix', '-R', regions_output_name, file_x], stdout=PIPE)
    #f = f.communicate()[0]
    #Remove universal_newlines=True if using python 2.7
    f = Popen(["/usr/local/bin/bcftools","query","-R", regions_output_name, file_x, "-f","%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%DS\t]\n"], stdout=PIPE, universal_newlines=True)
    f = f.communicate()[0]
    #If the tabix command finds one of the risk variants in that chunk:
    if f != '':
        #This will return the effect allele and effect for each line if variant is in the risk score and not a comment line, and then the dosage line (the := operator is specific to python 3.8, can replace fields with line.rstrip().split() which should work with any python, but will be slower as it does splitting multiple times
        test_results = np.asarray([flatten((score_file_dict[(fields[0] + ":" + fields[1] + ":" + fields[3] + ":" + fields[4])], fields[0:])) for line in f.rstrip().split("\n") if (fields := line.rstrip().split())[0][0] != "#" if (fields[0] + ":" + fields[1] + ":" + fields[3] + ":" + fields[4]) in score_file_dict])
        #Format of test_results is: effect allele, effect, chr, pos, variant_id, ref, alt, dosage*n_samples
        #G 0.2341 22 16050075 rs587697622 A G 0 0 0 0....
        if test_results.size != 0:
            #Where effect allele from risk score formula matches alternative allele, multiply directly
            matching_effect_allele = test_results[np.where(test_results[:,0] == test_results[:,6])]
            #[:, np.newaxis] this is needed to do the multipication element wise (first column * all dosages in row)
            matching_effect_allele = matching_effect_allele[:,1].astype(float)[:, np.newaxis] * matching_effect_allele[:,7:].astype(float)
            #Where effect allele matches reference allele, take 2-dosage, then multiply (so flip dosage to be for alternative allele)
            matching_reference_allele = test_results[np.where(test_results[:,0] == test_results[:,5])]
            matching_reference_allele = matching_reference_allele[:,1].astype(float)[:, np.newaxis] * (2 - matching_reference_allele[:,7:].astype(float))
            #Sum down columns
            dosage_scores_sum = np.sum(matching_reference_allele, axis=0) + np.sum(matching_effect_allele, axis=0)
            for x in range(len(dosage_scores_sum)):
                sample_score_dict[sample_id[x]] = sample_score_dict[sample_id[x]] + dosage_scores_sum[x]

print(timeit.default_timer() - start_time)

with open(args.output_file, 'w') as out:
    out.write("%s\t%s\n" % ("individual", "score"))
    for x in range(len(sample_id)):
        out.write("%s\t%.8f\n" % (sample_id[x], sample_score_dict[sample_id[x]]))
