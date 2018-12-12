

import sys, os, subprocess, gzip

from random import seed, sample
seed(42)

from genome_window_iter import genome_window_iter


# def basepair_gen(ref, smpl):
#     """
#     Generator for segments to write to each pseudohaploid
#     """
#     homozygote_chars = frozenset(list('ATGCNatgcn-'))
#     heterozygote_chars = frozenset(list('RYSWKMryswkm'))
       
#     iupac = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
#              'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT',
#              'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT',
#              'H': 'ACT', 'V': 'ACG', 'N': 'ATGC' }
#     # turn into tuples
#     for key, val in iupac.items():
#         iupac[key] = tuple(val)

#     Q = ord('Q')       
#     for r, s in zip(ref, smpl):
#         if s == Q:
#             yield (r, r)
#         elif chr(s) in homozygote_chars:
#             yield (s, s)
#         else:
#             assert chr(s) in heterozygote_chars, s
#             yield tuple(ord(u) for u in sample(iupac[chr(s)], 2))

def write_pseudo_haploid_samples(sample_file_name, ref_file_name, output_file_name1, output_file_name2):

    homozygote_chars = frozenset(list('ATGCNatgcn-'))
    heterozygote_chars = frozenset(list('RYSWKMryswkm'))
       
    iupac = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
             'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT',
             'K': 'GT', 'M': 'AC'}#, 
#             'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG'} #, 'N': 'ATGC' }
    # turn into tuples
    for key, val in iupac.items():
        iupac[key] = tuple(val)

    Q = ord('Q')   

    outfile_1 = gzip.open(output_file_name1, 'wb')
    outfile_2 = gzip.open(output_file_name2, 'wb')

    prev_smpl_chrom = None
    for smpl_window, ref_window in genome_window_iter(sample_file_name, ref_file_name, window_size=100000):

        smpl_chrom, smpl_start, smpl_end, smpl_seq = smpl_window
        ref_chrom, ref_start, ref_end, ref_seq = ref_window

        assert ref_chrom == smpl_chrom
        assert ref_start == smpl_start
        assert ref_end == smpl_end

        if not prev_smpl_chrom or smpl_chrom != prev_smpl_chrom:
            # fasta header
            header_template = prev_smpl_chrom and '\n>{}\n' or '>{}\n'
            outfile_1.write(header_template.format(smpl_chrom).encode())
            outfile_2.write(header_template.format(smpl_chrom).encode()) 

        prev_smpl_chrom = smpl_chrom

        out1 = list()
        out2 = list()
        for r, s in zip(ref_seq, smpl_seq):
            if s == 'Q':
                out1.append(r)
                out2.append(r)
            elif s in homozygote_chars:
                out1.append(s)
                out2.append(s)
            else:
                assert s in heterozygote_chars, s
                het = sample(iupac[s], 2)
                out1.append(het[0])
                out2.append(het[1])

        outfile_1.write(''.join(out1).encode())
        outfile_2.write(''.join(out2).encode())
            
    outfile_1.write(b'\n')
    outfile_2.write(b'\n')
            
# def write_pseudo_haploid_samples(sample_file_name, ref_file_name, output_file_name1, output_file_name2):
    
#     # load reference 
#     cmd = ref_file_name.endswith('.rz') and 'zcat' or 'cat'
#     p = subprocess.Popen([cmd, ref_file_name], 
#                          stdout=subprocess.PIPE, 
#                          stderr=subprocess.PIPE)
#     stdout, stderr = p.communicate()
#     if stderr: print(stderr)
#     smpl = dict()
#     for s in stdout.decode().split('>')[1:]:
#         header, s = s.split('\n', 1)
#         chrom = header.split()[0]
#         ref[chrom] = s.replace('\n', '').encode()
        
#     outfile_1 = gzip.open(output_file_name1, 'wb')
#     outfile_2 = gzip.open(output_file_name2, 'wb')

#     # chromosome names as strings
#     chromosome_names = [str(x) for x in range(1,23)] + ['X']


#     # generte pseudohaploid chromosomes
#     for chrom in chromosome_names:

#         # fasta header
#         outfile_1.write('>{}\n'.format(chrom).encode())
#         outfile_2.write('>{}\n'.format(chrom).encode())

#         buf1, buf2 = list(), list()
#         for i, (b1, b2) in enumerate(basepair_gen(ref[chrom], smpl[chrom])):

#             # collect 1000 segments before writing
#             buf1.append(b1)
#             buf2.append(b2)
            
#             # write segments in buffer
#             if not i % 10000:
#                 outfile_1.write(bytes(buf1))
#                 outfile_2.write(bytes(buf2))
# #                 outfile_1.writelines(buf1)
# #                 outfile_2.writelines(buf2)
#                 buf1, buf2 = [], []

#         # flush buffer
#         if len(buf1):
#             outfile_1.write(bytes(buf1))
#             outfile_2.write(bytes(buf2))
# #             outfile_1.writelines(buf1)
# #             outfile_2.writelines(buf2)

#         # trailing newline
#         outfile_1.write(b'\n')
#         outfile_2.write(b'\n')
  

if __name__ == "__main__":

    write_pseudo_haploid_samples(*sys.argv[1:])




