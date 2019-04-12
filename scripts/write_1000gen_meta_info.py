import sys, os

_, sample_name_file, sample_info_file, outdir = sys.argv

with open(sample_name_file) as f:
    vcf_samples = set(f.read().split())


n = 0
out_files = dict()
desc = dict()

all_males = open(os.path.join(outdir, 'all_males.txt'), 'w')
all_females = open(os.path.join(outdir, 'all_females.txt'), 'w')

with open(sample_info_file) as f:
    next(f) # skip header
    for line in f:
        sample, family_id, pop, pop_desc, gender, *rest = line.split(';')
        if sample not in vcf_samples:
            continue        
        desc[pop] = pop_desc

        if gender == 'male':
            print(sample, file=all_males)
        elif gender == 'female':
            print(sample, file=all_females)
        else:
            assert 0


        if (pop, gender) not in out_files:
            out_files[(pop, gender)] = open(os.path.join(outdir, '{}_{}.txt'.format(pop, gender)), 'w')
        print(sample, file=out_files[(pop, gender)])
        n += 1

with open(os.path.join(outdir, 'pop_names.tsv'), 'w') as f:
    for tup in desc.items():
        print(*tup, sep='\t', file=f)

assert n == len(vcf_samples)
