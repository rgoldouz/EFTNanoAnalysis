import subprocess
import sys
import os

infiles = sys.argv[11:]
print infiles
val = sys.argv[1:11]
print val
text = ''
text += '    TChain* ch    = new TChain("IIHEAnalysis") ;\n'
text += '    TChain ch1("meta") ;\n'
for fn in infiles:
    a, b = fn.split(':')
    text += '    ch->Add("' +  b + '");\n'
    text += '    ch1.Add("' +  b + '");\n'
text += '    ch ->AddFriend("meta");\n'
text += '    MyAnalysis t1(ch);\n'
text += '    t1.Loop("' + val[0]+'", "' + val[1] + '" , "'+ val[2] + '" , "'+ val[3] + '" , "'+ val[4] + '" , ' + val[5] + ' , '+ val[6] + ' , '+ val[7] + ' , '+ val[8] + ' , '+ val[9] +');\n'
SHNAME1 = 'main.C'
SHFILE1='#include "MyAnalysis.h"\n' +\
'int main(){\n' +\
text +\
'}'


open(SHNAME1, 'wt').write(SHFILE1)
#os.system("echo '" + SHFILE1+ "' > main.C && ls -l && cat main.C && root  -b -q -l main.so main.C && ls -l")
#os.system('echo ' + SHFILE1+ ' > main.C && ls -l && cat main.C && root -l main.so main.C && ls -l')
#os.system('cat main.C')
os.system('root -b -q -l libEFTGenReaderEFTHelperUtilities.so main.so main.C')
