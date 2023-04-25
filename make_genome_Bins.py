
import sys

genome_ver = sys.argv[1]
resolution = int(sys.argv[2])

chrsize_dict = {}    
f = open('./FAI/{}.fa.fai'.format(genome_ver), 'r')
for line in f:
	line = line.rstrip()
	linedata = line.split()
	chrsize_dict[linedata[0]] = int(linedata[1])
f.close()

for chrname in chrsize_dict.keys():
	o=open('./genome_Bin/{}.{}.{}.bin'.format(genome_ver, chrname, resolution),'w')
	i=0 
	size=int(chrsize_dict[chrname])
	while i<size:
		o.write('{}\t{}\t{}\n'.format(chrname, i, i+resolution))
		i+=resolution
	o.close()
#

