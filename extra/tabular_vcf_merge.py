def merge( split_files, output_file ):
	"""
	Multiple SAM files may each have headers. Since the headers should all be the same, remove
	the headers from files 1-n, keeping them in the first file only
	"""
	cmd = 'mv %s %s' % ( split_files[0], output_file )
	result = os.system(cmd)
	if result != 0:
		raise Exception('Result %s from %s' % (result, cmd))
	if len(split_files) > 1:
		cmd = 'egrep -hv "^#" %s >> %s' % ( ' '.join(split_files[1:]), output_file )
		result = os.system(cmd)
	if result != 0:
		raise Exception('Result %s from %s' % (result, cmd))
merge = staticmethod(merge)
