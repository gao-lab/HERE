import numpy as np
import sys
from collections import Counter

## all coordinates are 1-based, from 5' to 3' on reads
read_end_boundary_5_prime = 5
read_end_boundary_3_prime = int(sys.argv[1]) ## 90-10+1=81 for GSE44183, and (100-15)-10+1=76 for GSE36552

if [81, 76].__contains__(read_end_boundary_3_prime) == False:
    raise ValueError("Invalid read_end_boundary_3_prime: " + str(read_end_boundary_3_prime) + "; it should be 90-10+1=81 for GSE44183, and (100-15)-10+1=76 for GSE36552")

## output order: chr, pos, ref, total cov from mpileup, total nonN cov, ref+ cov, ref- cov, ref total cov, variant base, variant+ cov, variant- cov, variant total cov, variant read-end cov, variant read-middle cov
## for the following temp_line, we need ot set read_end_boundary_3_prime as 81 (this is a sample from GSE44183)
## read_end_boundary_3_prime = 81
## temp_line="chr1\t89923\ta\t92\t..T.T.T,T,T.,,,,T,,.,T,T..TTTtTTT.TT.tTTT,TtTT,.tTTTttttt.T.,,tt.tctTTt,,,tt.,,,T.t,.,Ttt,,.\t>=????@<?>??>?>>??>A>@???@???B??@@???A???@?@??=?@???A@AAA??A?>AA?A7A?=@???AA@?@?@?@?=><@@=><\t87,87,84,84,83,83,82,9,81,10,80,79,12,12,15,16,74,17,18,71,20,70,21,68,67,65,63,62,61,32,58,56,56,55,55,55,55,37,51,51,48,43,46,45,42,41,56,33,58,32,32,32,60,60,60,63,63,25,25,25,66,67,67,68,22,69,69,70,18,18,73,74,74,75,75,76,14,77,77,81,9,7,84,84,5,86,4,87,87,87,87,2"

mismatch_pairs_list = [['A', 'a'], ['C', 'c'], ['G', 'g'], ['T', 't']]

for temp_line in sys.stdin:
    ##print(temp_line.strip())
    temp_fields_list = temp_line.split(sep="\t")
    temp_counter = Counter(temp_fields_list[4])
    ## remove N's
    _ = temp_counter.pop('<', 0)
    _ = temp_counter.pop('>', 0)
    temp_total_nonN_count = sum(temp_counter.values())
    temp_ref_plusstrand_count = temp_counter.pop('.', 0)
    temp_ref_minusstrand_count = temp_counter.pop(',', 0)
    temp_ref_count = temp_ref_plusstrand_count + temp_ref_minusstrand_count
    temp_first_eight_fields_and_seps = "\t".join(temp_fields_list[0:4]) + "\t" + str(temp_total_nonN_count) + "\t" + str(temp_ref_plusstrand_count) + "\t" + str(temp_ref_minusstrand_count) + "\t" + str(temp_ref_count) + "\t"
    for temp_mismatch_base_plusstrand, temp_mismatch_base_minusstrand in mismatch_pairs_list:
        temp_variant_plusstrand_count = temp_counter.pop(temp_mismatch_base_plusstrand, 0)
        temp_variant_minusstrand_count = temp_counter.pop(temp_mismatch_base_minusstrand, 0)
        temp_variant_count = temp_variant_plusstrand_count + temp_variant_minusstrand_count
        if temp_variant_count > 0:
            temp_variant_readend_count = 0
            temp_variant_readmiddle_count = 0
            for temp_read_index, (temp_read_base_lower, temp_read_position) in enumerate(zip(temp_fields_list[4].lower(), np.fromstring(temp_fields_list[6], dtype="int", sep=","))):
                if temp_read_base_lower == temp_mismatch_base_minusstrand:
                    if (temp_read_position <= read_end_boundary_5_prime) or (temp_read_position >= read_end_boundary_3_prime):
                        temp_variant_readend_count = temp_variant_readend_count + 1
                    else:
                        temp_variant_readmiddle_count = temp_variant_readmiddle_count + 1
                ## print([temp_read_index, (temp_read_base_lower, temp_read_position, temp_variant_readend_count, temp_variant_readmiddle_count)])
            print(temp_first_eight_fields_and_seps +
                  temp_mismatch_base_plusstrand + "\t" +
                  str(temp_variant_plusstrand_count) + "\t" +
                  str(temp_variant_minusstrand_count) + "\t" +
                  str(temp_variant_count) + "\t" +
                  str(temp_variant_readend_count) + "\t" +
                  str(temp_variant_readmiddle_count) + "\t" 
                  )
