import shlex
import os
import sys
import csv
import statistics
from subprocess import Popen, PIPE
from joblib import Parallel, delayed
import multiprocessing
import re
import operator
import time
import random

_alg = "legacy parameter, not used"

# Produces the following output: UNIQUE_ID graph_name |V| |E| zoltan_cutn zoltan_stddev ratio_with_alg_dist ratio_with_zoltan

def generate_inps(parts, imbal, algname):
	coarsening_methods = ["agg", "rs_alg_dist"]
	mtx_files_all = os.listdir(".")
	mtx_files = list(filter(lambda x: ".mtx" in x, mtx_files_all))
	for f in mtx_files:
		fname = f.split(".")[0]
		for p in parts:
			for im in imbal:
				for alg in coarsening_methods:
					outfile_name = ".".join(["zdrive.inp", fname, alg, str(p), str(int(im*100))])
					with open(outfile_name, 'w') as outfile:
						print("Decomposition Method = hypergraph", file=outfile)
						print("Zoltan Parameters = HYPERGRAPH_PACKAGE=phg", file=outfile)
						print("Zoltan Parameters = lb_approach=partition", file=outfile)
						print("File Type = matrixmarket", file=outfile)
						print("File Name = ", fname, file=outfile)
						print("Parallel Disk Info = number=0", file=outfile)
						print("Zoltan Parameters = NUM_GLOBAL_PARTITIONS = ", p, file=outfile)
						print("Zoltan Parameters = PHG_COARSENING_METHOD=", alg, file=outfile)
						print("Zoltan Parameters = PHG_CUT_OBJECTIVE=HYPEREDGES", file=outfile)
						print("Zoltan Parameters = IMBALANCE_TOL=",im, file=outfile)
						print("ZOLTAN Parameter = SEED=",random.randint(1000,10000000), file=outfile)

parts = [2,4,8,16,32,64,128]

imbal = [1.03, 1.05, 1.1]
random.seed()
if len(sys.argv) == 2:
		if (sys.argv[1] != "agg" and sys.argv[1] != "rs_alg_dist"):
				print("Incorrect algorithm name. Only acceps \"agg\" and \"rs_alg_dist\"")
				sys.exit(-1)
		_alg = sys.argv[1]

generate_inps(parts, imbal, _alg);

# Get input files
