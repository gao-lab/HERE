import sys
import pandas as pd
miR_in_target_species_file = sys.argv[1]
All_file=sys.argv[2]
miR_in_all_species_file = sys.argv[3]

All = pd.read_csv(All_file, sep="\t")
miR_in_target_species = pd.read_csv(miR_in_target_species_file, sep="\t", header=None)
for i in range(miR_in_target_species.shape[0]):
    miR_family = miR_in_target_species.iloc[i, 0]
    tmp = All.loc[All['miR family'] == miR_family, 'Species ID']
    tmp = list(tmp.unique())
    tmp = [str(x) for x in tmp]
    miR_in_target_species.iloc[i, 2] = ";".join(tmp)
miR_in_target_species.to_csv(miR_in_all_species_file, index=False, header=False, sep="\t")
