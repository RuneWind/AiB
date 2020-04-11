import os

normal_files = os.listdir("trees/normal")
permuted_files = os.listdir("trees/permuted")

print(normal_files)
print(permuted_files)	


#for i in range(len(files)):
#	print("\n")
#	for j in range(i, i+1):
#		os.system("python rfdist.py trees/normal/" + files[i] + " trees/normal/" + files[j])