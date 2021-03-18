from gwf import Workflow

gwf = Workflow(
    defaults={'account': 'baboons'}
    )

gwf.target(name='vcfmerge',
     inputs=['results/analyzed_individuals.csv'], 
     outputs=['steps/vcf/nonafr_analyzed_individuals_chr7.vcf.gz'], 
     walltime='03:00:00', 
     memory='8g') << f"""
mkdir -p steps/vcf

tail -n +2  ~/simons/faststorage/people/kmt/results/analyzed_individuals.csv | grep -v Africa | cut -f 1 -d ',' | grep -f - ~/simons/faststorage/data/metadata/nature18964-s2-fixed-genders.csv | cut -f 3 -d ';' | awk '$1="/home/kmt/simons/faststorage/data/vcfs/"$1".annotated.nh2.variants.vcf.gz"' > steps/vcf/nonafr_analyzed_individuals_vcf_files.txt

bcftools merge --regions 7 --missing-to-ref -Oz -o steps/vcf/nonafr_analyzed_individuals_chr7.vcf.gz --file-list steps/vcf/nonafr_analyzed_individuals_vcf_files.txt

tabix steps/vcf/nonafr_analyzed_individuals_chr7.vcf.gz
"""

dedicated_indiv_list =[
'LP6005443-DNA_A06', # S_Greek-2, WestEurasia
'LP6005443-DNA_D06', # S_Korean-1, EastAsia
'LP6005519-DNA_D05', # S_Irula-2, SouthAsia
'LP6005443-DNA_D02', # S_Yakut-2, CentralAsiaSiberia
'LP6005443-DNA_F08', # S_Papuan-9, Oceania
'LP6005441-DNA_E10', # S_Pima-1, America
]
max_missing = 50000
smc_file_name_list = []
for dedicated_indiv in dedicated_indiv_list:

    smc_file_name = f'steps/smcpp/{dedicated_indiv}_{max_missing}.txt'
    smc_file_name_list.append(smc_file_name)

    gwf.target(name=f'vcf2smc_{dedicated_indiv.replace("-", "_")}',
        inputs=['steps/vcf/nonafr_analyzed_individuals_chr7.vcf.gz'], 
        outputs=[smc_file_name], 
        # walltime='3-00:00:00', 
        # memory='8g'
        ) << f"""
    mkdir -p steps/smcpp

#    SAMPLES=SS6004476,SS6004467,SS6004469,SS6004477,SS6004472,LP6007068-DNA_A01,SS6004468,SS6004474,LP6005519-DNA_D01,LP6005441-DNA_G06,LP6005443-DNA_G11,LP6005441-DNA_E10,LP6005677-DNA_E01,LP6005443-DNA_A12,LP6005677-DNA_D01,LP6005443-DNA_A03,LP6005442-DNA_F02,LP6005443-DNA_C03,LP6005443-DNA_D03,LP6005443-DNA_B03,LP6005443-DNA_C04,LP6005677-DNA_B02,LP6005441-DNA_E08,LP6005443-DNA_E05,LP6005443-DNA_D02,LP6005442-DNA_C07,LP6005443-DNA_G05,LP6005442-DNA_E07,LP6005519-DNA_B06,LP6005519-DNA_A06,LP6005441-DNA_G03,LP6005443-DNA_B01,LP6005441-DNA_C05,LP6005441-DNA_G05,LP6005441-DNA_C06,LP6005592-DNA_C02,LP6005442-DNA_C11,LP6005443-DNA_D06,LP6005443-DNA_E01,LP6005441-DNA_C08,LP6005443-DNA_E09,LP6005441-DNA_E09,LP6005443-DNA_F01,LP6005443-DNA_A07,LP6005443-DNA_H01,LP6005443-DNA_A02,LP6005443-DNA_B02,LP6005442-DNA_D01,LP6005443-DNA_C02,LP6005442-DNA_G01,LP6005592-DNA_H03,LP6005519-DNA_D06,LP6005592-DNA_B02,LP6005443-DNA_B08,LP6005443-DNA_C07,LP6005443-DNA_G07,LP6005443-DNA_H07,LP6005441-DNA_A10,LP6005443-DNA_E07,LP6005443-DNA_D07,LP6005443-DNA_C08,LP6005443-DNA_F07,LP6005443-DNA_A08,LP6005443-DNA_F08,LP6005441-DNA_D01,LP6005441-DNA_C01,LP6005442-DNA_G09,LP6005519-DNA_G03,LP6005519-DNA_H03,LP6005441-DNA_D03,LP6005441-DNA_C03,LP6005441-DNA_E03,LP6005441-DNA_F05,LP6005441-DNA_E05,LP6005519-DNA_C05,LP6005519-DNA_D05,LP6005441-DNA_E06,LP6005519-DNA_A04,LP6005519-DNA_B04,LP6005519-DNA_E05,LP6005443-DNA_C09,LP6005443-DNA_D09,LP6005519-DNA_G04,LP6005519-DNA_H04,LP6005441-DNA_C07,LP6005519-DNA_E04,LP6005519-DNA_F04,LP6005441-DNA_C10,LP6005592-DNA_B04,LP6005442-DNA_A12,LP6005519-DNA_A05,LP6005519-DNA_B05,LP6005441-DNA_G11,LP6005519-DNA_C04,LP6005519-DNA_D04,LP6005442-DNA_C02,LP6005442-DNA_D02,LP6005441-DNA_A01,LP6005442-DNA_G02,LP6005519-DNA_F03,LP6005441-DNA_C02,LP6005441-DNA_E02,LP6005441-DNA_A06,LP6005442-DNA_A03,LP6005442-DNA_B03,LP6005442-DNA_D03,LP6005443-DNA_H05,LP6005441-DNA_G04,LP6005442-DNA_E10,LP6005442-DNA_H03,LP6005442-DNA_G03,LP6005592-DNA_A02,LP6005441-DNA_A05,LP6005442-DNA_B04,LP6005442-DNA_A04,LP6005442-DNA_G07,LP6005443-DNA_A06,LP6005442-DNA_A08,LP6005442-DNA_C04,LP6005443-DNA_B10,LP6005592-DNA_E01,LP6005592-DNA_G03,LP6005442-DNA_F04,LP6005442-DNA_E04,LP6005442-DNA_G04,LP6005443-DNA_E10,LP6005443-DNA_F10,LP6005441-DNA_C09,LP6005592-DNA_B03,LP6005592-DNA_E02,LP6005441-DNA_G10,LP6005592-DNA_D01,LP6005592-DNA_D04,LP6005441-DNA_C11,LP6005442-DNA_A11,LP6005519-DNA_C03,LP6005519-DNA_D03,LP6005677-DNA_C03,LP6005441-DNA_G12,LP6005592-DNA_G01

    SAMPLES=`tail -n +2  ~/simons/faststorage/people/kmt/results/analyzed_individuals.csv | grep -v Africa | cut -f 1 -d ',' | grep -f - ~/simons/faststorage/data/metadata/nature18964-s2-fixed-genders.csv | cut -f 3 -d ';' | tr '\n' ',' | sed 's/.$//'`

    singularity run terhorst/smcpp:latest vcf2smc --cores 4 --missing-cutoff {max_missing} -d {dedicated_indiv} {dedicated_indiv} steps/vcf/nonafr_analyzed_individuals_chr7.vcf.gz {smc_file_name} 7 nonAfr:$SAMPLES
    """

gwf.target(name='estimate',
     inputs=smc_file_name_list, 
     outputs=['model.final.json'], 
     walltime='4-00:00:00', 
     memory='16g'
     ) << f"""

singularity run docker://terhorst/smcpp:latest estimate -o steps/smcpp --cores 4 --timepoints 35 4e4 --spline pchip 1.247e-08 {' '.join(smc_file_name_list)}

"""

gwf.target(name='plot',
     inputs=['model.final.json'], 
     outputs=['nonAfr_joint_demography.png'], 
    #  walltime='4-00:00:00', 
    #  memory='16g'
     ) << f"""

docker run --rm -v $PWD:/mnt terhorst/smcpp:latest plot -g 29 -c nonAfr_joint_demography.png model.final.json
"""




# import os, sys
# cores=4,
#      walltime='4-00:00:00', 
#      memory='16g') << f"""

# singularity run docker://terhorst/smcpp:latest estimate -o steps/smcpp --cores 4 --timepoints 35 4e4 --spline pchip 1.247e-08 {' '.join(smc_file_name_list)}

# """

# gwf.target(name='plot',
#      inputs=['model.final.json'], 
#      outputs=['nonAfr_joint_demography.png'], 
#     #  walltime='4-00:00:00', 
#     #  memory='16g'
#      ) << f"""

# docker run --rm -v $PWD:/mnt terhorst/smcpp:latest plot -g 29 -c nonAfr_joint_demography.png model.final.json
# """



# demography_chr7 = scipy.interpolate.PchipInterpolator(generations_chr7, Ne_chr7, extrapolate=True)
# demography_chrX = scipy.interpolate.PchipInterpolator(generations_chrX, Ne_chrX, extrapolate=True)

# plt.plot(generations_chr7, [demography_chrX(g) / demography_chr7(g) for g in generations_chr7])

