
import os, sys, tempfile, shutil

_, truncation_sample = sys.argv
truncation_sample = int(truncation_sample)

for file_path in sys.stdin:
	file_path = file_path.strip()
	# get info about orig time stamp
	stats = os.stat(file_path)

	# make a named tmp file
	tmpfile_handle, tmpfile_name = tempfile.mkstemp(text=True)

	# write the relevant part of the file into the tmp file
	tmpfile = os.fdopen(tmpfile_handle, 'w')
	with open(file_path) as f:
		for line in f:
			stage, iteration, *_ = line.split()
			if stage == 'resample' and int(iteration) > truncation_sample:
				break
			tmpfile.write(line)
	tmpfile.close()

	#  move file to overwrite other
	shutil.move(tmpfile_name, file_path)

	# change time stamp back to what is was:
	os.utime(file_path, ns=(stats.st_atime_ns, stats.st_mtime_ns))

