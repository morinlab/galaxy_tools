import argparse

if __name__ == '__main__':

	desc = "Convert Strelka Makefile into a sub-Makefile used STRICTLY for parallelization by chromosome in Galaxy"
	parser = argparse.ArgumentParser(description=desc)
	parser.add_argument(
		"-m", "--makefile",
		required=True,
		help="Input Must be a Strelka Makefile"
		)
	parser.add_argument(
		"-c", "--chrom",
		type=str,
		required=True,
		help="Chromosome present in the Strelka Makefile"
		)
	parser.add_argument(
		'-o', '--output',
		required=True,
		help="Output Makefile"
		)
	args = parser.parse_args()

	makefile_in = open(args.makefile, 'r')
	makefile_out = open(args.output, 'w')

	count = 0
	for line in makefile_in:
		if count > 0:
			makefile_out.write(line)
			count = count - 1
			if count == 0:
				makefile_out.write('\n')
		elif count < 0:
			count = count + 1
			continue
		elif line.startswith('script_dir'):
			makefile_out.write(line)
		elif line.startswith('call_script'):
			makefile_out.write(line)		
		elif line.startswith('filter_script'):
			makefile_out.write(line)
		elif line.startswith('finish_script'):
			makefile_out.write(line)
		elif line.startswith('script_dir'):
			makefile_out.write(line)
		elif line.startswith('config_file'):
			makefile_out.write(line)
		elif line.startswith('analysis_dir'):
			makefile_out.write(line)
		elif line.startswith('results_dir'):
			makefile_out.write(line)		
		elif line.startswith('get_chrom_dir'):
			makefile_out.write(line)
		elif line.startswith('get_chrom_task'):
			makefile_out.write(line)
		elif line.startswith('get_bin_task'):
			makefile_out.write(line)
			makefile_out.write('\n')
		elif line.startswith('all:'):
			makefile_out.write(line)
			makefile_out.write('\n')
		elif line.startswith('$(finish_task):'):
			makefile_out.write(line)
			makefile_out.write('\techo "DONE"\n')
			makefile_out.write('\n')
		elif line.startswith('_'.join(['chrom',args.chrom,'bin'])):
			makefile_out.write(line)
			count = 3	
		elif line.startswith('_'.join(['chrom',args.chrom,'task'])):
			makefile_out.write(line)
			count = 3
		elif line.startswith('chrom'):
			count = -3

	makefile_in.close()
	makefile_out.close()