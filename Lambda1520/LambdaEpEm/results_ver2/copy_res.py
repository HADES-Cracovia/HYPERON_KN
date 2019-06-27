

directory="../output_ver2/" #oryginal data directory
fname="../files_list_k.dat" #file with list of expected files
import subprocess
print("Copy and hadd files from \"results_ver2\" directory")
import time
import os

with open(fname) as f:
    content = f.readlines()
content = [x.strip() for x in content]

for k in content:   #take every name from vector content
    subdir= ''.join([i for i in k if i.isdigit()])
    opendir=directory+'0'+subdir+"/"
    
    #run hadd for all files from directory
    bashCommand = "hadd -f "+k+" "+ opendir +"*.root"
    print(bashCommand)
    os.system(bashCommand)

print("coped following files:")
for k in content:
    print(k)
