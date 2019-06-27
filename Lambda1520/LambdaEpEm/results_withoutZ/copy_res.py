directory="/lustre/nyx/hades/user/iciepal/Lambda1520_ic/" #oryginal data directory
fname="../files_list_k.dat" #file with list of expected files
import subprocess
print("Copy and hadd files from Iza directory")
import time
import os

with open(fname) as f:
    content = f.readlines()
content = [x.strip() for x in content]

for k in content:   #take every name from vector content
    subdir=result = ''.join([i for i in k if i.isdigit()])
    opendir=directory+"channel"+subdir+"/"
    
    #run hadd for all files from directory
    bashCommand = "hadd -f "+k+" "+ opendir +"*.root.out"
    print(bashCommand)
    os.system(bashCommand)

print("coped following files:")
for k in content:
    print(k)
