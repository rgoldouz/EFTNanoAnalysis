import os, sys, argparse, uproot
import subprocess

def GetListOfWCs(fname):
  ''' Retruns a list of the WC names from WCnames, (retruns [] if not an EFT sample) '''
  wc_names_lst = []
  tree = uproot.open(f'{fname}:Events')
  if 'WCnames' not in tree.keys():
    wc_names_lst = []
  else:
    wc_info = tree['WCnames'].array(entry_stop=1)[0]
    for idx,i in enumerate(wc_info):
      h = hex(i)[2:]                                 # Get rid of the first two characters
      wc_fragment = bytes.fromhex(h).decode('utf-8') # From: https://stackoverflow.com/questions/3283984/decode-hex-string-in-python-3
      # The WC names that are longer than 4 letters are too long to be encoded in a 64-bit integer:
      #   - They're instead stored in two subsequent entries in the list
      #   - This means that the decoded names in wc_info go like this [... 'ctlT' , '-i' ...]
      #   - The leading '-' indicates the given fragment is the trailing end of the previous WC name
      #   - The following logic is supposed to put those fragments back together into the WC name
      if not wc_fragment.startswith("-"):
        wc_names_lst.append(wc_fragment)
      else:
        leftover = wc_fragment[1:]                    # This leftover part of the WC goes with the previous one (but get rid of leading '-')
        wc_names_lst[-1] = wc_names_lst[-1]+leftover  # So append this trailing fragment to the leading framgenet to reconstruct the WC name
  return wc_names_lst

import topcoffea.modules.fileReader as fr
fpath = "/afs/crc.nd.edu/user/r/rgoldouz/MakeLobsterJobs/UL/mgprod/lobster_workflow/ul_cfgs/NAOD-00000.root"
print(fr.GetListOfWCs(fpath))
