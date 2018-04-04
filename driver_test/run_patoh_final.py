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

# Produces the following output: UNIQUE_ID graph_name |V| |E| zoltan_cutn

def generate_inps(parts, imbal):
	mtx_files_all = os.listdir(".")
	mtx_files = list(filter(lambda x: ".mtx" in x, mtx_files_all))
	for f in mtx_files:
		fname = f.split(".")[0]
		alg = "patoh"
		for p in parts:
			for im in imbal:
				outfile_name = ".".join(["zdrive.inp", fname, alg, str(p), str(int(im*100))])
				with open(outfile_name, 'w') as outfile:
						print("Decomposition Method = hypergraph", file=outfile)
						print("Zoltan Parameters = HYPERGRAPH_PACKAGE=patoh", file=outfile)
						print("Zoltan Parameters = lb_approach=partition", file=outfile)
						print("File Type = matrixmarket", file=outfile)
						print("File Name = ", fname, file=outfile)
						print("Parallel Disk Info = number=0", file=outfile)
						print("Zoltan Parameters = NUM_GLOBAL_PARTITIONS = ", p, file=outfile)
						print("Zoltan Parameters = IMBALANCE_TOL=",im, file=outfile)

def get_first_word_after_substring(output, str_to_look_for):
	outsplit = output.split(str_to_look_for, 1)
	if len(outsplit) < 2:
		return None
	return str((re.match(b'\w+', outsplit[1])).group(0))

def get_first_int_after_substring(output, str_to_look_for):
	outsplit = output.split(str_to_look_for)
	if len(outsplit) < 2:
		return None
	outclean = []
	for s in outsplit:
		outclean.append(s.strip())

	return int(float((re.match(b'\d+', outclean[-1])).group(0)))

def get_first_num_after_substring(output, str_to_look_for):
	outsplit = output.split(str_to_look_for)
	if len(outsplit) < 2:
		return None
	outclean = []
	for s in outsplit:
		outclean.append(s.strip())

	return float((re.match(b'\d+.\d+', outclean[-1])).group(0))

def run_for_file(filename, num_iter):

	cutns = []
	alg_name = ""
	num_parts = -1	

	for i in range(num_iter):
		cmd = "./zdrive.exe " + filename

		# Run zdrive
		process = Popen(shlex.split(cmd), stdout=PIPE)
		(output, err) = process.communicate()
		exit_code = process.wait()

		#Process output
		alg_name = get_first_word_after_substring(output, b'HYPERGRAPH_PACKAGE = ')
		num_parts = get_first_int_after_substring(output, b'Zoltan_LB_Eval_HG  Part count:')
		imbalance = get_first_num_after_substring(output, b'IMBALANCE_TOL[0] =')
		cutn = get_first_num_after_substring(output, b'CUTN (Sum_edges( (#parts(edge)>1)*ewgt )):')
		nvtx = get_first_int_after_substring(output, b'RS_VERTEX_NUM=')
		nedge = get_first_int_after_substring(output, b'RS_HEDGE_NUM=')
		npins = get_first_int_after_substring(output, b'RS_PINS_NUM=')

		if (alg_name == None) or (num_parts == None) or (cutn == None) or (imbalance == None) or (cutn == 0):
			print("Something went wrong")
			continue
		cutns.append(cutn)

	graphname = (filename.split("."))[-4]
	unique_id = str(graphname) + str(num_parts) + str(imbalance)

	if len(cutns) == 0:
		result = [unique_id, graphname, nvtx, nedge, npins, num_parts, imbalance, -1]
	else:
		result = [unique_id, graphname, nvtx, nedge, npins, num_parts, imbalance, min(cutns)]

	tmp_out_filename = "output_for_inp_" + filename + ".csv"
	with open(tmp_out_filename, 'w') as tmp_out:
	    tmp_w = csv.writer(tmp_out, delimiter=';')
	    tmp_w.writerow(result)
	return result 
	

parts = [2,4,8,16,32,64,128]

outfilename = 'output_'
if len(sys.argv) >= 2:
	outfilename += sys.argv[1]

outfilename += '.csv'

if len(sys.argv) <= 2:
	imbal = [1.03, 1.05, 1.1]
else:
	try:
		imbal = [float(sys.argv[2])]
	except ValueError:
		print("Incorrect imbalance encountered, exiting")
		sys.exit(-1)

generate_inps(parts, imbal);


# Get input files

inp_files_all = os.listdir(".")
inp_files = list(filter(lambda x: "zdrive.inp." in x, inp_files_all))

num_cores = 8 
results_unsorted = Parallel(n_jobs=num_cores)(delayed(run_for_file)(f, 10) for f in inp_files)
results = sorted(results_unsorted, key=operator.itemgetter(1, 2, 3))

with open(outfilename, 'w') as csvfile:
	out = csv.writer(csvfile, delimiter=';')
	out.writerow(["patoh", "patoh", "patoh", "patoh", "patoh", "patoh", "patoh", "patoh"])
	for r in results:
		out.writerow(r)

