'''
python tau.py  -f ../data/test.out
'''

import pyalps
from report import report

import argparse
parser = argparse.ArgumentParser(description='')
parser.add_argument("-fileheaders", nargs='+', default="params", help="fileheaders")
args = parser.parse_args()

resultFiles = []
for fileheader in args.fileheaders:
    resultFiles += pyalps.getResultFiles(prefix=fileheader)

resultFiles = list(set(resultFiles))

for filename in resultFiles:

    obslist = ['PertOrder','M2', 'IntE','IntE2', 'Energy']

    print filename 
    report(filename, obslist)
    print "#######################################################################"
