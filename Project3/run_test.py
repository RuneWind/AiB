import os
from Bio import SeqIO
import subprocess

d = sorted(os.listdir("./testseqs"))

exact_score = []
approx_score = []

for file in d:
	print(file)
	result_e = subprocess.check_output("python sp_exact_3.py sub_m.txt 5 ./testseqs/" + file, shell=True, text=True)
	exact_score.append(int(result_e.strip()))
	print("Exact: ", exact_score)
	#print(result_e.strip())

	result_a = subprocess.check_output("python sp_approx.py sub_m.txt 5 ./testseqs/" + file + " True", shell=True, text=True)
	approx_score.append(int(result_a.strip()))
	print("Approx: ", approx_score)
	#print(result_a.strip())

print(exact_score)
print(approx_score)