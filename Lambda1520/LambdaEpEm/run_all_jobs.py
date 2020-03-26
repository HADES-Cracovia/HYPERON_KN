directory="/lustre/nyx/hades/user/iciepal/Lambda1520_ic/" #oryginal data directory
fname="list_allfiles_newGeometryM3" #file with list of expected files
import subprocess
print("run job for all channels\n")
import time
import os

with open(fname) as f:
    content = f.readlines()
content = [x.strip() for x in content]

print("File list: \n")
print(content)
#os.system("scancel -u knowakow")
#print("scancel -u knowakow")
while(1 != int(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,))):
        #print('waiting to finish list: {}'.format(k))
        print('still {} jobs to the end'.format(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,)))
        time.sleep(30*4)


for k in content:   #take every name from vector content
    print('try to run files from list:\n{}'.format(k))
    while(1 != int(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,))):
        print('waiting to finish list: {}'.format(k))
        print('still {} jobs to the end'.format(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,)))
        time.sleep(30)
    bashCommand = "./run_job.py "+k+" -d output_newGeometryM3_ver3_pi0"
    print(bashCommand)
    os.system(bashCommand)

#print("coped following files:")
#for k in content:
 #   print(k)
