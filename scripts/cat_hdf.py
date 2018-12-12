
import pandas
import sys
file_name = sys.argv[1]

#print(pandas.read_hdf(file_name))
for row in pandas.read_hdf(file_name).itertuples():
    print(*row, sep='\t')
