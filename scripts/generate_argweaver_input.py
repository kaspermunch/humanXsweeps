#!/usr/bin/env python

import sys
import os
from optparse import OptionParser

from genome_window_iter import genome_window_iter
import simons_meta_data


def file_base_name(p):
    b = os.path.basename(p)
    name, suf = os.path.splitext(b)
    while suf:        
        name, suf = os.path.splitext(name)
    return name

def main():

    usage="""
%prog [options] [inputFileGlob [outputFile]]

"""
    parser = OptionParser(usage=usage, version="%prog 1.0")

    parser.add_option("-w", "--windowsize",
                      dest="windowsize",
                      type="int",
                      default=None,
                      help="window size to write to different files")
    parser.add_option("-o", "--outputdir",
                      dest="outputdir",
                      type="str",
                      default=None,
                      help="output dir")


    (options, args) = parser.parse_args()

    file_name_list = args

    assert options.windowsize is not None
    assert options.outputdir is not None

    if not os.path.exists(options.outputdir):
        os.makedirs(options.outputdir)

#    full_file_name_list = [l.strip() for l in file_name_list]

    meta_data, populations, regions = simons_meta_data.get_meta_data()

#    window_iter = genome_window_iter(*file_name_list, window_size=options.windowsize, chunk_size=options.chunksize)
    window_iter = genome_window_iter(*file_name_list, window_size=options.windowsize)

    for window in window_iter:

#        names, starts, ends, _ = list(zip(*window))
        names, starts, ends, seqs = list(zip(*window))

        assert names[1:] == names[:-1]
        assert starts[1:] == starts[:-1]
        assert ends[1:] == ends[:-1]

        outfile = os.path.join(options.outputdir, "{}-{:09d}-{:09d}.fa".format(names[0], starts[0], ends[0]))

        with open(outfile, 'w') as f:
            for (name, start, end, seq), file_name in zip(window, file_name_list):
                print(">{}\n{}\n".format(file_base_name(file_name), seq), file=f)


if __name__ == "__main__":
    main()



