import numpy as np
import sys

## temp_line="chr8\t28190741\tN\t24\tAGaAaAaAAAAaaGAAAaAnGaaa\tGC@GG:<GGGGGGGGFGGG#GGGG\t115,111,86,85,85,73,73,67,61,60,57,55,48,47,39,38,36,28,21,19,16,50,10,10"

reference_and_N_chars_list = [".", ",", ">", "<"]

for temp_line in sys.stdin:
    temp_fields_list = temp_line.split(sep="\t")
    temp_tens_list = []
    temp_ones_list = []
    temp_base_quality_passed_list = []
    temp_base_quality_literal = temp_fields_list[5].strip()
    for temp_i, temp_base_quality in enumerate(temp_base_quality_literal):
        temp_i_corrected = temp_i + 1
        temp_tens = int(temp_i_corrected/10)
        temp_ones = int(temp_i_corrected) - temp_tens * 10
        temp_tens_list.append(str(temp_tens))
        temp_ones_list.append(str(temp_ones))
        temp_base_quality_passed_list.append(str(int(ord(temp_base_quality) >= 58)))
        
    print("".join(temp_tens_list))
    print("".join(temp_ones_list))
    print(temp_fields_list[4])
    print(temp_base_quality_literal)
    print("".join(temp_base_quality_passed_list))        
