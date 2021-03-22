from gwf import Workflow

gwf = Workflow(
    defaults={'account': 'baboons'}
    )

###################################################

dedicated_indiv_list =[
'LP6005443-DNA_A06', # S_Greek-2, WestEurasia
'LP6005443-DNA_D06', # S_Korean-1, EastAsia
'LP6005519-DNA_D05', # S_Irula-2, SouthAsia
'LP6005443-DNA_D02', # S_Yakut-2, CentralAsiaSiberia
# 'LP6005443-DNA_F08', # S_Papuan-9, Oceania
# 'LP6005441-DNA_E10', # S_Pima-1, America
]
max_missing = 100000
mutation_rate = 1.45e-08 # 5e-10 * 29
#mutation_rate = 1.247e-08
generation_time = 29

###################################################

# run vcfmerge
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

# run vcf2sms
smc_file_name_list = []
for dedicated_indiv in dedicated_indiv_list:

    smc_file_name = f'steps/smcpp/vcf2smc/{dedicated_indiv}_{max_missing}.txt'
    smc_file_name_list.append(smc_file_name)

    gwf.target(name=f'vcf2smc_{dedicated_indiv.replace("-", "_")}',
        inputs=['steps/vcf/nonafr_analyzed_individuals_chr7.vcf.gz'], 
        outputs=[smc_file_name], 
        walltime='02:00:00', 
        memory='8g'
        ) << f"""
    mkdir -p steps/smcpp/vcf2smc

    SAMPLES=`tail -n +2  ~/simons/faststorage/people/kmt/results/analyzed_individuals.csv | grep -v Africa | cut -f 1 -d ',' | grep -f - ~/simons/faststorage/data/metadata/nature18964-s2-fixed-genders.csv | cut -f 3 -d ';' | tr '\n' ',' | sed 's/.$//'`

    singularity run docker://terhorst/smcpp:latest vcf2smc --cores 4 --missing-cutoff {max_missing} -d {dedicated_indiv} {dedicated_indiv} steps/vcf/nonafr_analyzed_individuals_chr7.vcf.gz {smc_file_name} 7 nonAfr:$SAMPLES
    """
    
# run estimate and plot   
prefix = f"pchip_{max_missing}_{'_'.join(dedicated_indiv_list)}"    
gwf.target(name='estimate_pchip111',
     inputs=smc_file_name_list, 
     outputs=[f'steps/smcpp/{prefix}/model.final.json', f'steps/smcpp/{prefix}/{prefix}.png'], 
     walltime='4-00:00:00', 
     cores=10,
     memory='16g'
     ) << f"""
mkdir -p steps/smcpp/{prefix}
singularity run docker://terhorst/smcpp:latest estimate -o steps/smcpp/{prefix} --cores 10 --timepoints 35 4e4 --spline pchip {mutation_rate} {' '.join(smc_file_name_list)}
singularity run docker://terhorst/smcpp:latest plot -g {generation_time} -c steps/smcpp/{prefix}/{prefix}.png steps/smcpp/model.final.json
"""

# run estimate and plot   
prefix = f"step_{max_missing}_{'_'.join(dedicated_indiv_list)}"    
gwf.target(name='estimate_step111',
     inputs=smc_file_name_list, 
     outputs=[f'steps/smcpp/{prefix}/model.final.json', f'steps/smcpp/{prefix}/{prefix}.png'], 
     walltime='4-00:00:00', 
     cores=10,
     memory='16g'
     ) << f"""
mkdir -p steps/smcpp/{prefix}
singularity run docker://terhorst/smcpp:latest estimate -o steps/smcpp/{prefix} --cores 10 --timepoints 35 4e4 1.247e-08 {' '.join(smc_file_name_list)}
singularity run docker://terhorst/smcpp:latest plot -g {generation_time} -c steps/smcpp/{prefix}/{prefix}.png steps/smcpp/{prefix}/model.final.json
"""



# demography_chr7 = scipy.interpolate.PchipInterpolator(generations_chr7, Ne_chr7, extrapolate=True)
# demography_chrX = scipy.interpolate.PchipInterpolator(generations_chrX, Ne_chrX, extrapolate=True)

# plt.plot(generations_chr7, [demography_chrX(g) / demography_chr7(g) for g in generations_chr7])

