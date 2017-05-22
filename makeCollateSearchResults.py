#!/usr/bin/python
#
###########################################
#
# File: makeCollateSearchResults.py
# Author: Ra Inta
# Description:
# Created: October 19, 2016
# Last Modified: October 21, 2016
#
###########################################


# The C versions of Element Tree are faster (but not always supported)
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
try:
    from xml.etree.cElementTree import Element, SubElement, Comment, tostring
except ImportError:
    from xml.etree.ElementTree import Element, SubElement, Comment, tostring

import bz2
from intervaltree import Interval,IntervalTree
import os
from sys import argv
import re
from math import ceil, floor
import glob

from CasACommon import CondorSubFile,assert_no_files_matching,indent


############################################################
# 0) Read in some basic setup information
############################################################

JOBS_PER_SUBDIR=250  # For the moment this is hard-coded here (following the style of two other scripts in the pipeline).

# If we have an input argument, then take it as the number of parallel jobs to 
# cut the collate job into, otherwise take it as 200
N_PARALLEL_JOBS = int(argv[1]) if len(argv)>1 else 200 


SCRIPTS = os.getcwd()


# Get parameters from search_setup.xml
with open('../../search_setup.xml','r') as searchSetup:
    setup_tree = ET.parse( searchSetup )
    setup_root = setup_tree.getroot()

for search_params in setup_root.iter('search'):
    search_band = float( search_params.find('band').text )
    start_freq = float( search_params.find('freq').text )
    Tspan = float( search_params.find('span_time').text )

for ul_params in setup_root.iter('upper_limit'):
    ul_band = float( ul_params.find('band').text )

NUM_UL_BANDS = int( ceil( search_band / ul_band ) )

############################################################
# 1) Take in search results
############################################################

with bz2.BZ2File('../../search_bands.xml.bz2', 'rb') as searchBandsFile:
    search_tree = ET.parse( searchBandsFile )
    search_root = search_tree.getroot()

## Walk through the XML tree for job, 2F and f0; don't pick up the vetoed bands.
search_segments = []

for jobNumber in search_root.iter('search_band'):
    search_jobNum = int( jobNumber.find('job').text )
    search_startFreq = float( jobNumber.find('freq').text )
    search_endFreq = search_startFreq + float( jobNumber.find('band').text )
    search_segments.append( (search_startFreq, search_endFreq, search_jobNum) )

# Take the search_segments for this job as a cue
NUM_SEARCH_BANDS = search_segments[-1][-1] + 1


# check current directory is empty of previous job files
assert_no_files_matching(['condor.*', 'search.sub.*'])

## TODO explicitly set the subCollate search job numbers
# here, and make the subsequent script take in corresponding arguments?

indent(search_root)

def writeSearchBandJobIdx(jobIdx, jobStartIdx, jobEndIdx):
    """docstring for writeSearchBandJobIdx"""
    # Write search_bands.xml
    search_sub_file = "search_bands_"+str(jobIdx)+".xml"
    search_bands = ET.Element("search_bands")
    for eachEl in search_root[jobStartIdx:jobEndIdx+1]:
        search_bands.append( eachEl )
    search_bands_xml = ET.ElementTree( search_bands )
    search_bands_xml.write(search_sub_file, xml_declaration=True, encoding='UTF-8', method='html')
    os.system("bzip2 -f " + search_sub_file)

## Write upper limits XML
#indent(upper_limit_root)
#upper_limit_bands_xml = ET.ElementTree(upper_limit_root)
#upper_limit_bands_xml.write("upper_limit_bands_py.xml", xml_declaration=True, encoding='UTF-8', method='html' )

# Get list of how jobs should be cut up. This can be a little tricky as we don't
# want overlapping search jobs, so it's not evenly spaced.
jobStep = int(ceil(1.0*NUM_SEARCH_BANDS/(N_PARALLEL_JOBS)))
jobsList = range( 0, NUM_SEARCH_BANDS, jobStep)
jobsList.append( NUM_SEARCH_BANDS )



# Cut full collate job into N_PARALLEL_JOBS jobs
# TODO balance load by examining filesizes of search_results.txt?
for jobIdx in range(N_PARALLEL_JOBS):
    CondorSubObj = CondorSubFile('collateSearchResults.py')
    #CondorSubObj.args.append( jobIdx )
    CondorSubObj.args = [jobIdx]
    #CondorSubObj.args.append( N_PARALLEL_JOBS )
    #CondorSubObj.args.append( NUM_SEARCH_BANDS )
    CondorSubObj.fileName = 'collate_search_results.sub.' + str(jobIdx)
    CondorSubObj.output = 'condor.out.' + str(jobIdx)
    CondorSubObj.error = 'condor.err.' + str(jobIdx)
    CondorSubObj.write()
    # Cut up search_bands.xml.bz2 to search_bands_jobIdx.xml.bz2
    job_start = jobsList[jobIdx]
    job_end = jobsList[jobIdx+1] - 1
    writeSearchBandJobIdx(jobIdx, job_start, job_end)
    #print(str(job_start) + "   " + str(job_end))


