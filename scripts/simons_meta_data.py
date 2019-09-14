import os
from collections import defaultdict
import numpy
from pathlib import Path

###########################################################
# various meta data settings of my own that apply across analyses:
###########################################################

# group1_pure_africans = ['S_Ju_hoan_North-1', 'S_Ju_hoan_North-2', 'S_Ju_hoan_North-3', 
#                         'S_Khomani_San-1', 'S_Khomani_San-2', 
#                         'S_Mbuti-1', 'S_Mbuti-2', 'S_Mbuti-3']

# group2_pure_africans = ['S_BantuHerero-1', 'S_BantuHerero-2', 'S_BantuKenya-1', 
#                         'S_BantuKenya-2', 'S_BantuTswana-1', 'S_BantuTswana-2', 
#                         'S_Biaka-1', 'S_Biaka-2', 
#                         'Dinka-1', 'Dinka-2', 
#                         'S_Esan-1', 'S_Esan-2', 
#                         'S_Gambian-1', 'S_Gambian-2', 
#                         'S_Igbo-1', 'S_Igbo-2', 
#                         'S_Kongo-2', 
#                         'S_Lemande-1', 
#                         'S_Lemande-2', 
#                         'S_Luhya-1', 'S_Luhya-2', 
#                         'S_Luo-1', 'S_Luo-2', 
#                         'S_Mandenka-1', 'S_Mandenka-2', 
#                         'S_Mende-1', 'S_Mende-2', 
#                         'S_Mozabite-1', 'S_Mozabite-2', 
#                         'S_Saharawi-1', 'S_Saharawi-2', 
#                         'S_Yoruba-1', 'S_Yoruba-2']

# pure_africans = group1_pure_africans + group2_pure_africans

# Excluded individuals due to bad quality (these are taken out of the "individuals" meta data):
excluded_individuals = ['S_Palestinian-2', # elise and moi found these had strange coverage patterns
                        'S_Naxi-2',        # elise and moi found these had strange coverage patterns
                        'S_Lezgin-1',      # Genetic sex is not assigned
                        "S_Finnish-1", "S_Finnish-2", "S_Mansi-1", "S_Mansi-2", "S_Palestinian-2"]
                        # the last five are also excluded in Mallick et at.: "5 samples based on missing X chromosome 
                        # data in an initial processing, for themselves or a second sample"

# African populations with backflow
# (these are removed in the subsequent notebook analysis and in computing argweaver pruned_tmrca_half):
#excluded_populations = ['Masai', 'Somali' ,'Mozabite', 'Saharawi', 'Luhya', 'Luo', 'BantuKenya', 'Gambian']

# exclude sub-Saharan African populations:
excluded_populations = ['Masai', 'Somali', # show some non-African component in Structure plot in main paper (Laurits does that too)
                        'Mozabite', 'Saharawi', # known non-african ancestry and Neanderthal admxiture
                        'Luhya', 'Luo', 'BantuKenya', # 1000 genomes study sugggest that Kenyans (Luhya) have some backflow. (Laurits info)
                        'Gambian'] # Not sure what the argument is here, but Laurits excludes it too...

###########################################################
# read in meta data to get sex of each individual
###########################################################

def get_meta_data(meta_data_dir=None, include_ust_ishim=False):

    if meta_data_dir is not None:
        d = meta_data_dir
    else:
        d = Path('/home/kmt/simons/faststorage/data/metadata')

    meta_data_file_name = d / 'nature18964-s2.csv'

    meta_data = dict()
    with open(str(meta_data_file_name), 'r') as f:
        keys = f.readline().split(';')
        for l in f:
            d = dict(zip(keys, l.split(';')))
            sample_id = d['Sample ID (SGDP)']
            sample_pop = d['Population ID']
            if sample_id not in excluded_individuals and sample_pop not in excluded_populations \
                    and d['Embargo level (X=Fully Public, Y=Signed Letter)'] == 'X':             
                try:
                    d['Longitude'] = float(d['Longitude'].replace(',', '.'))
                    d['Latitude'] =  float(d['Latitude'].replace(',', '.'))
                except ValueError:
                    d['Longitude'] = numpy.nan
                    d['Latitude'] =  numpy.nan
                meta_data[sample_id] = d

    if include_ust_ishim:
        d = dict()
        d['Longitude'] = 71.167194
        d['Latitude'] = 57.692596
        d['Population ID'] = 'Ust_Ishim'
        d['Region'] = 'Ust_Ishim'
        d['Genetic sex assignment'] = 'XY'
        meta_data['Ust_Ishim'] = d

    # dicts of samples from each population and region
    populations = defaultdict(list)
    regions = defaultdict(list)
    for key in meta_data:
        pop = meta_data[key]['Population ID']
        if pop not in excluded_populations:
            populations[pop].append(key)

            reg = meta_data[key]['Region']
            regions[reg].append(key)

    return meta_data, populations, regions

if __name__ == "__main__":

    meta_data, populations, regions = get_meta_data()

    print("Regions:")
    for r in regions:
        print('\t', r)

    print("Populations:")
    for p in sorted(populations):
        print('\t', p)

    print("Individuals:")
    for i in sorted(meta_data):
        print('\t', i)
