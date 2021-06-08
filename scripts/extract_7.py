#!/usr/bin/env python

import sys, gzip

input_file_name, output_file_name = sys.argv[1], sys.argv[2]

with gzip.open(input_file_name, 'rb') as input_file:
    for l in input_file:
        if l.decode().startswith('>7'):
            with open(output_file_name, 'wb') as output_file:
                output_file.write(l)
                output_file.write(input_file.readline())
            break
