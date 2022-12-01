import numpy as np
import sys

## temp_line="chr1\t77770\tT\t27\t..............,c.......,.,.\tSRR893046.9330128,SRR893046.9557551,SRR893046.5050182,SRR893046.10161386,SRR893046.1419897,SRR893046.3918386,SRR893046.6079033,SRR893046.4667129,SRR893046.12410010,SRR893046.1786604,SRR893046.371057,SRR893046.11857026,SRR893046.5596746,SRR893046.4544441,SRR893046.13485688,SRR893046.5697229,SRR893046.6796255,SRR893046.5029489,SRR893046.7685179,SRR893046.2815428,SRR893046.11288029,SRR893046.3484461,SRR893046.13317237,SRR893046.6450030,SRR893046.8060316,SRR893046.3417115,SRR893046.15320651"

reference_and_N_chars_list = [".", ",", ">", "<"]

for temp_line in sys.stdin:
    temp_fields_list = temp_line.split(sep="\t")
    temp_first_four_fields_and_seps = "\t".join(temp_fields_list[0:4]) + "\t"
    for temp_index, (temp_base, temp_read_name) in enumerate(zip(temp_fields_list[4], temp_fields_list[5].split(",")), start=1):
        if reference_and_N_chars_list.__contains__(temp_base) == False:
            print(temp_first_four_fields_and_seps + str(temp_index) + "\t" + temp_base + "\t" + temp_read_name)
