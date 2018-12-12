



import gzip, sys


# def _window_iter(file_name, window_size=None, chunk_size=1):

#     assert window_size

#     if file_name.endswith('.gz') or file_name.endswith('.rz'):
#         open_file = gzip.open
#     else:
#         open_file = open

#     with open_file(file_name, 'rb') as f:

#         block, data = '', ''
#         chrom, seq = None, None
#         idx = 0
#         while True:

#             # get data
#             if len(block) < 10 * window_size: # block quaratied to not to 
#                                               # accumulate with many small fasta entries

#                 # read a chunk of the file
#                 data = f.read(chunk_size).decode()

#                 # at end of last fasta
#                 if not data and len(block.replace('\n', '')) < window_size and block.find('>') < 0:
#                     if block.replace('\n', ''):
#                         seq = block.replace('\n', '')
#                         yield chrom, idx, idx+len(seq), seq
# #                        idx +=  len(seq)
#                     return

#                 # add data to buffer
#                 block += data

#             # see if there is the start of a header in the block
#             header_start = block.find('>')
#             print(block[header_start])
#             if header_start >= 0:

#                 # see if there is some sequence from the prev entry we need to flush
#                 # before parsing the found header
#                 if block[:header_start].replace('\n', ''):
#                     assert chrom is not None

#                     # write sequence blocks
#                     buf = block[:header_start].replace('\n', '')
#                     while buf:
#                         seq = buf[:window_size]
#                         yield chrom, idx, idx+len(seq), seq
#                         idx += len(seq)
#                         buf = buf[window_size:]

#                 # find the end of the header
#                 header_end = block.find('\n', header_start)
#                 # keep reading in chunks until we find it
#                 while header_end < 0:
#                     data = f.read(chunk_size).decode()
#                     block += data
#                     header_end = block.find('\n', header_start)
#                     assert data or header_end, "file truncated in the middle of a header"

#                 # slice out header
#                 header = block[header_start:header_end]

#                 # record current chrom
#                 chrom = header[1:].split()[0]
#                 # reset start index
#                 idx = 0

#                 # cut header of start of buffer
#                 block = block[header_end+1:]

#                 continue

#             # write sequence blocks
#             buf = block.replace('\n', '')
#             while len(buf) >= window_size:
#                 seq = buf[:window_size]
#                 yield chrom, idx, idx+len(seq), seq
#                 idx += len(seq)
#                 buf = buf[window_size:]
#             block = buf

def _reader(f, chunksize):

    buf = f.read(chunksize).decode()
    while buf:
        for b in buf:
            yield b
        buf = f.read(chunksize).decode()


def _window_iter(file_name, window_size=None, skipped=[]):

    skipped = set(skipped)

    assert window_size

    if file_name.endswith('.gz') or file_name.endswith('.rz'):
        open_file = gzip.open
    else:
        open_file = open

    with open_file(file_name, 'rb') as f:

        char = _reader(f, 100000)

        chrom, seq = None, []
        idx = 0
        max_header = 10000
        while True:

            #c = f.read(1).decode()
            c = next(char)

            # skip newlines
            if c == '\n':
                continue

            # flush at end of file
            if not c:
                if seq:
                    yield chrom, idx, idx + len(seq), ''.join(seq)
                return

            # start of header
            elif c == '>':
                # flush
                if seq:
                    yield chrom, idx, idx + len(seq), ''.join(seq)
                    seq = []

                # get header
                header = ''
                chrom = ''
                for i in range(max_header):
                    #h = f.read(1).decode()
                    h = next(char)
                    if h == '\n':
                        chrom = header.split()[0]
                        idx = 0
                        break
                    header += h
                assert chrom, 'header end not found'

            # else:
            elif chrom not in skipped:
                # add to sequence buffer
                seq.append(c)

                # write chunk of seq if necessary
                if len(seq) >= window_size:
                    yield chrom, idx, idx + window_size, ''.join(seq[:window_size])
                    seq = seq[window_size:]                    
                    idx += window_size



#def genome_window_iter(*genome_file_names, window_size=None, chunk_size=100000):
def genome_window_iter(*genome_file_names, window_size=None, skip=[]):
    assert len(genome_file_names)
    assert type(window_size) == int
#    assert type(chunk_size) == int
    if len(genome_file_names) == 1:
#        return _window_iter(genome_file_names[0], window_size=window_size, chunk_size=chunk_size)
        return _window_iter(genome_file_names[0], window_size=window_size, skipped=skip)
    else:
#        return zip(*(_window_iter(f, window_size=window_size, chunk_size=chunk_size) for f in genome_file_names))
        return zip(*(_window_iter(f, window_size=window_size, skipped=skip) for f in genome_file_names))


if __name__ == "__main__":

    for window in genome_window_iter(*sys.argv[1:], window_size=10):

        print(window)

        names, starts, ends, seqs = list(zip(*window))

        assert names[1:] == names[:-1]
        assert starts[1:] == starts[:-1]
        assert ends[1:] == ends[:-1]










