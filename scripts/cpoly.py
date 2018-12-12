

import simons_meta_data

individuals, populations, regions = simons_meta_data.get_meta_data()

f = open('samples.ind', 'w')

print('Chimp', file=f)

for indiv in individuals:
    chromotype = individuals[indiv]['Genetic sex assignment']
    sex = chromotype == 'XY' and 'M' or 'F'
    pop = individuals[indiv]['Population ID']
    print(indiv, sex, pop, file=f)

        
