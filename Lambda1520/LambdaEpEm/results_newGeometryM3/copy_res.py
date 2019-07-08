

directory="../output_newGeometryM3/" #oryginal data directory
fname="../files_list_M3.dat" #file with list of expected files
#fname="file5758_list.dat"

import subprocess
print("Copy and hadd files from \"output_new_geometry_2\" directory")
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
