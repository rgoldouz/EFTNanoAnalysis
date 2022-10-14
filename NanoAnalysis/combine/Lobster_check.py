import subprocess
import sys
import os

quark = ['DU','DC','SU','SC','BU','BC']
Couplings = ['cT','cS']
theoryXS =[[32.53,11.32],[2.79,0.98],[8.19,2.8],[0.7,0.26],[3.24,1.11],[0.28,0.10]]


val = sys.argv[1:4]
print val
Y = val[2].split('_')
cardNameB1=''
#if 'Mu' in val[1]:
#    for year in Y:
#        cardNameB1 += val[0] +'_'+val[1]+'_emu_'+year+'_' + 'llB1.txt '+ val[0] +'_'+val[1]+'_mumu_'+year+'_' + 'llB1.txt '
#else:
#    for year in Y:
#        cardNameB1 += val[0] +'_'+val[1]+'_ee_'+year+'_' + 'llB1.txt ' + val[0] +'_'+val[1]+'_emu_'+year+'_' + 'llB1.txt '

for year in Y:
    cardNameB1 += val[0] +'_'+val[1]+'_ee_'+year+'_' + 'llB1.txt ' + val[0] +'_'+val[1]+'_emu_'+year+'_' + 'llB1.txt '+ val[0] +'_'+val[1]+'_mumu_'+year+'_' + 'llB1.txt '

rmax= 1400.0 / (1000*theoryXS[quark.index(val[1][1:3])][Couplings.index(val[0])])
rmin = -1 * rmax


print cardNameB1
os.system('cp CombinedFilesBNV/* .')
os.system('combineCards.py ' + cardNameB1 + ' > ' +val[0] +'_'+val[1]+'_'+val[2]  + '_com.txt')
#os.system('combine --run blind -M  AsymptoticLimits '         +val[0] +'_'+val[1]+'_'+val[2]  + '_com.txt > ' + val[0] +'_'+val[1]+'_'+val[2] +'_results.tex' )
os.system('combine -M  AsymptoticLimits '          +val[0] +'_'+val[1]+'_'+val[2]  + '_com.txt --rMax ' + str(rmax) + ' > ' + val[0] +'_'+val[1]+'_'+val[2] +'_results.tex' )
os.system('text2workspace.py  '                   +val[0] +'_'+val[1]+'_'+val[2] + '_com.txt -m 125')
#expected impact plot
#os.system('combineTool.py -M Impacts -d '         +val[0] +'_'+val[1]+'_'+val[2] + '_com.root -m 125 --doInitialFit --robustFit 1 --rMin -3 --rMax 3 -t -1 --expectSignal=1')
#os.system('combineTool.py -M Impacts -d '         +val[0] +'_'+val[1]+'_'+val[2] + '_com.root -m 125 --robustFit 1 --doFits --parallel 4 --rMin -3 --rMax 3 -t -1 --expectSignal=1 --parallel 8')
#os.system('combineTool.py -M Impacts -d '         +val[0] +'_'+val[1]+'_'+val[2] + '_com.root -m 125 -o impacts.json --rMin -3 --rMax 3 -t -1 --expectSignal=1')
os.system('combineTool.py -M Impacts -d '         +val[0] +'_'+val[1]+'_'+val[2] + '_com.root -m 125 --doInitialFit --robustFit 1 --rMin -3 --rMax 3 -t -1 --expectSignal=0')
os.system('combineTool.py -M Impacts -d '         +val[0] +'_'+val[1]+'_'+val[2] + '_com.root -m 125 --robustFit 1 --doFits --parallel 4 --rMin -3 --rMax 3 -t -1 --expectSignal=0 --parallel 8')
os.system('combineTool.py -M Impacts -d '         +val[0] +'_'+val[1]+'_'+val[2] + '_com.root -m 125 -o impacts.json --rMin -3 --rMax 3 -t -1 --expectSignal=0')
if len(Y)==1:
    os.system('plotImpacts.py -i impacts.json -o impacts --max-pages 1 --label-size 0.03 --cms-label ,' + val[0] +'_'+val[1]+'_'+val[2])
else:
    os.system('plotImpacts.py -i impacts.json -o impacts --max-pages 1 --label-size 0.03 --cms-label ,' + val[0] +'_'+val[1]+'_2016-2018')
os.system('mv impacts.pdf ' + val[0] +'_'+val[1]+'_'+val[2] + '_Expected_impacts.pdf')

#Observed impact plot
os.system('combineTool.py -M Impacts -d '         +val[0] +'_'+val[1]+'_'+val[2] + '_com.root -m 125 --doInitialFit --robustFit 1  --rMin ' + str(rmin/3.0) + ' --rMax ' + str(rmax/3.0))
os.system('combineTool.py -M Impacts -d '         +val[0] +'_'+val[1]+'_'+val[2] + '_com.root -m 125 --robustFit 1 --doFits --parallel 4 --parallel 8  --rMin ' + str(rmin/3.0) + ' --rMax ' + str(rmax/3.0))
os.system('combineTool.py -M Impacts -d '         +val[0] +'_'+val[1]+'_'+val[2] + '_com.root -m 125 -o impacts.json  --rMin ' + str(rmin/3.0) + ' --rMax ' + str(rmax/3.0))
if len(Y)==1:
    os.system('plotImpacts.py -i impacts.json -o impacts --max-pages 1 --label-size 0.025 --cms-label ,' + val[0] +'_'+val[1]+'_'+val[2])
else:
    os.system('plotImpacts.py -i impacts.json -o impacts --max-pages 1 --label-size 0.025 --cms-label ,' + val[0] +'_'+val[1]+'_2016-2018')
os.system('mv impacts.pdf ' + val[0] +'_'+val[1]+'_'+val[2] + '_Observed_impacts.pdf')





