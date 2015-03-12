import argparse
import math
import pysam

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
    	help="Input BAM file"
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
    parser.add_argument(
        "--chromosome",
        action="store_true",
        required=False,
        help="Add the Chromosome to the Interval File"
        )
    parser.add_argument(
        "--start_position",
        action="store_true",
        required=False,
        help="Add the Start Position to the Interval File"
        )
    parser.add_argument(
        "--end_position",
        action="store_true",
        required=False,
        help="Add the End Position to the Interval File"
        )

    args = parser.parse_args()
    data = {}
    alignment = pysam.AlignmentFile(args.input, "rb")
    header = alignment.header
    chunk_size = args.chunk_size

    for seq in header['SQ']:
        rname = seq['SN']
        length = seq['LN']
        if not rname in header:
            data[rname] = []
        
        if args.mode == 'by_rname':
            chunk_size = int(length)

        for i in range(int(math.ceil(float(length) / float(chunk_size)))):
            if (i+1) * chunk_size > length:
                data[rname].append([rname,(i * chunk_size) + 1, length])
            else:
                data[rname].append([rname,(i * chunk_size) + 1, (i+1) * chunk_size])

    output = open(args.output,'w')

    for seq in header['SQ']:
        rname = seq['SN']
        for interval in data[rname]:
            line = ''
            if args.chromosome:
                line = ''.join([line,'%s\t' % interval[0]])
            if args.start_position:
                line = ''.join([line,'%s\t' % interval[1]])
            if args.end_position:
                line = ''.join([line,'%s\t' % interval[2]])
            output.write(''.join([line[:-1], '\n']))