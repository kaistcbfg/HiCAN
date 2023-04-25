
import numpy as np
import cooler

import pickle
import gzip
import sys

clrfile = sys.argv[1]
resolution = int(sys.argv[2])
chrlist_filepath = sys.argv[3]
outfilename = sys.argv[4]

c500k = cooler.Cooler("{}::/resolutions/{}".format(clrfile,resolution))

f = open(chrlist_filepath)
chrlist = [i.rstrip() for i in f.readlines()]
f.close()

rows = []
for idx1, i in enumerate(chrlist):

	tmp = []
	balance_type = False

	for idx2, j in enumerate(chrlist):

		if idx1 == idx2:
			mat = c500k.matrix(balance=balance_type).fetch(i)
		elif idx1 < idx2:
			mat = c500k.matrix(balance=balance_type).fetch(j,i)
		elif idx1 > idx2:
			mat = c500k.matrix(balance=balance_type).fetch(j,i)
		#
		tmp.append(mat)

	#
	rowmat = np.concatenate(tmp,axis=0)
	rows.append(rowmat)
#
 
#from iced import normalization
#totalmat = normalization.ICE_normalization(totalmat)
#totalmat = np.nan_to_num(totalmat, nan=0, posinf=0, neginf=0)
totalmat = np.concatenate(rows,axis=1)

with gzip.open(outfilename, 'wb') as f: pickle.dump(totalmat,f)

