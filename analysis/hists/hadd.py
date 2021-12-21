import sys
import os
import subprocess
import readline
import string
import glob
sys.path.append('/afs/crc.nd.edu/user/r/rgoldouz/BNV/analysis/bin')
import Files_BNV
SAMPLES = {}
SAMPLES.update(Files_BNV.mcBNV_samples)



os.system('rm *.root')
dist = "/hadoop/store/user/rgoldouz/FullProduction/BNVAnalysis/Analysis_" 

for key, value in SAMPLES.items():
    year = value[3]
    hadd='hadd ' + key + '.root '
    for filename in os.listdir(dist + key):
        hadd += dist + key + '/' + filename + ' '
    os.system(hadd)


#for key, value in SAMPLES.items():
#    os.system('cp ' + glob.glob(dist  + key + '/*.root')[0] + ' ' + key + '.root')

