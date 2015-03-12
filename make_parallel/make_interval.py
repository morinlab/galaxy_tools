import argparse
import math

if __name__ == "__main__":
    desc = "Create BAM interval file for easily parallelism in Galaxy"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
    	"--mode", "-m", 
    	required=True,
    	choices=['by_rname', 'by_chunk']
    	)
    parser.add_argument(
    	"--input", "-i",
    	type=str,
    	required=True,
    	help="Input FASTA file"
    	)
    parser.add_argument(
    	"--output", "-o", 
    	type=str,
    	required=True,
    	help="Output BED file"
    	)
    parser.add_argument(
    	"--chunk_size", "-s", 
    	type=int,
    	required=False,
    	default=50000000,
    	help="The number of bases to include in a chunk interval, required when 'by_chunks' specified"
    	)
    args = parser.parse_args()
    data = {}
    rnames = list()
    reference = open(args.input, 'r')
    output = open(args.output, 'w')
    current_rname = ''

    for line in reference:
        line = line[:-1]
        line = line.replace(' ', '')
        if line.startswith('>'):
            current_rname = line[1:]
            rnames.append(current_rname)
            print current_rname
            data[current_rname] = {}
            data[current_rname]['GENOME'] = 0
            data[current_rname]['KNOWN'] = 0
            continue
        else:
            data[current_rname]['GENOME'] += len(line)
            line = line.replace('N','')
            data[current_rname]['KNOWN'] += len(line)

    if args.mode == 'by_chunk':
    	for rname in rnames:
    		for i in range(int(math.ceil(float(data[rname]['GENOME']) / float(args.chunk_size)))):
    			if (i+1) * args.chunk_size > data[rname]['GENOME']:
    				line = "%s\t%s\t%s\t%s\t%s\n" % (rname, (i * args.chunk_size) + 1, data[rname]['GENOME'], data[rname]['KNOWN'], data[rname]['GENOME'])
    				output.write(line)
    			else:
    				line = "%s\t%s\t%s\t%s\t%s\n" % (rname, (i * args.chunk_size) + 1, (i+1) * args.chunk_size, data[rname]['KNOWN'], data[rname]['GENOME'])
    				output.write(line)

    else:
    	for rname in rnames:
    		line = "%s\t%s\t%s\t%s\t%s\n" % (rname, 1, data[rname]['GENOME'], data[rname]['KNOWN'], data[rname]['GENOME'])
    		output.write(line)

   	output.close()