import sys 
import os 
import subprocess 
import readline 
import string 

mcBNV_samples={'BNV_TT_TSCE': [['rgoldouz/FullProduction/BNV/ntuple_BNV_TT_TSCE'], 'mc', 'none', '2017', 'none', '1', '41.53', '200000.0', '1', '40.0'], 'BNV_ST_TDUE': [['rgoldouz/FullProduction/BNV/ntuple_BNV_ST_TDUE'], 'mc', 'none', '2017', 'none', '1', '41.53', '200000.0', '1', '40.0'], 'BNV_TT_TSUE': [['rgoldouz/FullProduction/BNV/ntuple_BNV_TT_TSUE'], 'mc', 'none', '2017', 'none', '1', '41.53', '200000.0', '1', '40.0'], 'BNV_ST_TSCE': [['rgoldouz/FullProduction/BNV/ntuple_BNV_ST_TSCE'], 'mc', 'none', '2017', 'none', '1', '41.53', '200000.0', '1', '40.0'], 'BNV_ST_TSUE': [['rgoldouz/FullProduction/BNV/ntuple_BNV_ST_TSUE'], 'mc', 'none', '2017', 'none', '1', '41.53', '200000.0', '1', '40.0'], 'BNV_TT_TDUE': [['rgoldouz/FullProduction/BNV/ntuple_BNV_TT_TDUE'], 'mc', 'none', '2017', 'none', '1', '41.53', '200000.0', '1', '40.0']}
