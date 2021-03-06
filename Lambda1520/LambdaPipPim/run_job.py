#!/usr/bin/env python3

#----------------------- AUTHOR: J. Berger-Chen, July 2014 -------------------------------
#CALL EXAMPLE:
#python run_dst.py /hera/hades/dst/jul14/online/187/root/ test/ la_187_test 100 La --nFiles=1
#python run_dst.py /hera/hades/dst/jul14/online/188/root/ k0_188/ k0_188 K0 10000000
#python run_dst.py /hera/hades/dst/jul14/online/187/root/ la_187/ la_187 La 10000000

# modified by R. Lalik

import os
import sys
import glob
#sys.path.append('/hera/hades/user/klapidus/online_rec/python/argparse-1.2.1')
import argparse
import shlex, subprocess

def_time = 15
def_mem = "500mb"
def_script='job_script.sh'

parser=argparse.ArgumentParser(description='Submit dst analysis to GSI batch farm')
parser.add_argument('arguments', help='list of arguments', type=str, nargs='+')
parser.add_argument('-p', '--part', help='partition', type=str, default="main")
parser.add_argument('-f', '--file', help='input is single file', action='store_true', default=False)
parser.add_argument('-e', '--events', help='number of events per file to be processed',type=int, default=1000000)
parser.add_argument('-t', '--time', help='time need to finish task', default=def_time, type=int)
parser.add_argument('-m', '--mem', help='requested memory', default=def_mem, type=str)
parser.add_argument('-s', '--script', help='execute script', default=def_script)
parser.add_argument('-d', '--directory', help='output directory', default='.', type=str)
parser.add_argument('-n', '--files',help='number of files you want to proceed',type=int,default=-1)
parser.add_argument('-l', '--split',help='number of files you want to proceed',type=int,default=0)
args=parser.parse_args()

print(args)

NFILES = args.files

submissiondir=os.getcwd()+'/'
tpl_resources='--time={0:1d}:{1:02d}:00 --mem-per-cpu={2:s} -p main'
jobscript=submissiondir+args.script

if __name__=="__main__":
    i=0

    resources = tpl_resources.format(int(args.time/60), args.time % 60, args.mem)
    lines = []

    deps = [None] * args.split
    for arg in args.arguments:

        if args.file is False:
            f = open(arg)
            lines = f.readlines()
            f.close()
        else:
            lines = [ arg ]

        jobid_list = ""

        for entry in lines:
            entry = entry.strip()

            i += 1

            logfile = submissiondir + '/log/slurm-%j.log'

            events = args.events

            #print('submitting file: ',entry)

            if os.path.isfile(logfile):
                os.remove(logfile)

            command = "-o {:s} {:s} -J {:s} --export=\"pattern={:s},events={:d},odir={:s}\"".format(logfile, resources, os.path.split(entry)[1], entry, events, args.directory)
            if (args.split > 0 and deps[i % args.split] != None):
                command += " -d afterany:"+deps[i % args.split]
            command += " " + jobscript

            print('sbatch ' + command)
            job_command='sbatch ' + command
            proc = subprocess.Popen(shlex.split(job_command), stdout=subprocess.PIPE, shell=False)
            (out, err) = proc.communicate()
            if (out[0:19] == b'Submitted batch job'):
                jobid = out[20:-1].decode()
                if (args.split > 0):
                    deps[i % args.split] = jobid
            else:
                print("Job failed with error:")
                print(err)

        print(i, ' entries submitted')
   
