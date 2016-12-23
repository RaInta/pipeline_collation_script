#!/usr/bin/python
#
###########################################
#
# File: testFunctions.py
# Author: Ra Inta
# Description:
# Created: October 20, 2016
# Last Modified: 20161212, R.I. 
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
#from intervaltree import Interval,IntervalTree
import os
import sys
#import re
from math import ceil, cos
import glob


SCRIPTS = os.path.dirname(sys.argv[0])        


##################################################
### The following will eventually be       #######
### refactored---perhaps a separate text file  ###
### TODO this is ugly because they're globals! ###
##################################################

#LIGO_USER = "ra.inta"
#ACCT_TAG = "ligo.prod.o1.cw.directedisolatedother.coherent"
LIGO_USER = "ben.owen"
ACCT_TAG = "ligo.prod.o1.cw.directedisolated.coherent"

##################################################

############################################################
###   Check no files match the fileName  ###################
############################################################
def assert_no_files_matching( fileName ):
    """Checks that there are no files matching the (globbed list) input"""
    #TODO refactor this as try...catch
    for fileGlobIdx in fileName:
        for fileIdx in glob.glob( fileGlobIdx ):
            if os.path.isfile( fileIdx ):
                print(os.getcwd() + " is not empty of files matching " + fileGlobIdx )
                break
############################################################

############################################################
###   Create a Condor submit file as a class    ############
############################################################

class CondorSubFile:
    def __init__(self, executable, fileName='sngl_job.sub', args=[], queue=1, output='condor.out', error='condor.err'):
        """
        This class produces objects capable of producing Condor submit files.

        Sensible defaults are used; if you only want to create a single job, the
        only required input parameter is the name of the executable, e.g.:

            submitJob = CondorSubFile('Search.py')
            submitJob.fileName = 'search.sub'

        The optional args parameter is a (space delimited) list. This is ignored
        if not explicitly set. Otherwise, the list is appended to the submit file.
        e.g.:

            submitJob.args = ['--penguinType','Emperor','--curvature','11.3']

        The only method is a write() call to create the submit files after setting
        the appropriate parameters, such as the fileName. i.e.:

            submitJob.write()

        Will write out the submit file described.

        The idea is that a simple
        loop can be used to tweak parameters (such as job number) easily.
        """
        self.executable = executable
        self.fileName = fileName
        self.args = args
        self.queue = queue
        self.output = output
        self.error = error
        # TODO Let's get rid of these globals --- maybe put in search_setup.xml
        global LIGO_USER
        global ACCT_TAG

    def write(self):
        """The only method for this class; writes out the Condor submit file"""
        # The boilerplate stuff
        header_str = "# global options set by " +__file__+ "\n\n"
        header_str += "universe = vanilla\n"
        header_str += "initialdir = " + os.getcwd() + "\n"
        header_str += "getenv = true\n"
        header_str += "notification = never\n"
        header_str += "log = condor.log\n"
        header_str += "request_cpus = 1\n"
        header_str += "request_memory = 4 GB\n"
        header_str += "request_disk = 4 GB\n"
        header_str += "accounting_group = " + ACCT_TAG + "\n"
        header_str += "accounting_group_user = " + LIGO_USER + "\n"
        header_str += "\n"
        # Job specific options
        body_str = "# this job's options and commands\n\n"
        body_str += "executable = " + os.path.join(SCRIPTS, self.executable) + "\n"
        if len(self.args):
            body_str += "arguments = "
            body_str += " ".join([str(x) for x in self.args] )
            body_str += "\n"
        body_str += "output = " + self.output + "\n"
        body_str += "error = " + self.error + "\n"
        body_str += "queue " + str(self.queue) + "\n"
        with open(self.fileName,'w') as subFileName:
            subFileName.write( header_str + body_str )
############################################################



############################################################
###    Read setup XMLs into an object       ################
############################################################
class SetupXML:
    def __init__(self, fileName='search_setup'):
        """docstring for ReadSetupXML
        For example:
            searchSetup = SetupXML('search_setup')
            searchSetup.search.band = 200
            searchSetup.upper_limit.band = 1.0
            searchSetup.write()
        will write out the search_setup.xml, populated with the appropriate
        values.
            searchSetup.read()
            """
        # Get parameters from search_setup.xml
        self.fileName = fileName
        def read(self):
            """docstring for read"""
            pass
        with open('search_setup.xml','r') as searchSetup:
            setup_tree = ET.parse( searchSetup )
            setup_root = setup_tree.getroot()
        self.search
        def write(self):
            """docstring for write"""
            pass
        pass
############################################################

############################################################
# Test for possiblity of BZ2 zipped search results
############################################################
def openFile(fileName):
    if fileName[-4:] == '.bz2':
        return bz2.BZ2File(fileName, "r")
    else:
        return open(fileName, "r")
############################################################

#for search_params in setup_root.iter('search'):
#    search_band = float( search_params.find('band').text )
#    start_freq = float( search_params.find('freq').text )
#    Tspan = float( search_params.find('span_time').text )
#
#for ul_params in setup_root.iter('upper_limit'):
#    ul_band = float( ul_params.find('band').text )
#
#
#
#
#def read_setup_xmls():
#    """docstring for read_setup_xmls"""
#    pass
#
# TODO Port the following to Python:
## read setup XML files
#sub read_setup_xmls {
#
#    my %toread = @_;
#
#    my %setup;
#
#    # read setup XML files
#    foreach my $xmlfile ((SETUP_XML, SFT_STRETCH_XML, COMP_COST_XML)) {
#	next if !(-f $xmlfile);
#
#	# read XML file
#	my $twig = XML::Twig->new();
#	$twig->safe_parsefile($xmlfile) or
#	    croak "XML::Twig->safe_parsefile failed: $@";
#
#	# iterate over first-level elements
#	foreach my $elt1 ($twig->root->children()) {
#
#	    # iterate over second-level elements
#	    foreach my $elt2 ($elt1->children()) {
#
#		# parse element into hash
#		$setup{$elt1->name}->{$elt2->name} = $elt2->text;
#
#	    }
#
#	}
#
#    }
#
#    # iterate over requested values
#    foreach my $name (keys %toread) {
#
#	# get keys
#	my ($key1, $key2) = split /:/, $name;
#
#	# croak if keys does not exist
#	croak "Couldn't find '$name' in setup XML files"
#	    if !defined($setup{$key1}->{$key2});
#
#	# set referenced variable to value
#	${$toread{$name}} = $setup{$key1}->{$key2};
#
#    }
#
#}

# TODO Port the following to Python:
## write a setup XML file
#sub write_setup_xml {
#
#    my ($xmlfile, %towrite) = @_;
#
#    my %xml;
#
#    # iterate over requested values
#    foreach my $name (keys %towrite) {
#
#	# get keys
#	my ($key1, $key2) = split /:/, $name;
#
#	# create key and value in XML structure
#	$xml{$key1}->{$key2} = $towrite{$name};
#
#    }
#
#    # write XML setup file
#    my $xs = XML::Simple->new(
#			      XMLDecl  => 1,
#			      RootName => 'setup',
#			      NoAttr   => 1
#			      );
#    open FILE, ">$xmlfile" or croak "Couldn't open '$xmlfile': $!";
#    print FILE $xs->XMLout(\%xml);
#    close FILE;
#
#}

##################################################





############################################################
# ElementTree doesn't pretty-print XML! Need to add indents.
# Stolen directly from:
# http://stackoverflow.com/questions/749796/pretty-printing-xml-in-python
############################################################


def indent(elem, level=0):
    """This is a function to make the XML output pretty, with the right level
    of indentation"""
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

############################################################


############################################################
# Get template covering band
############################################################
def template_covering_band(loudestLine, Tspan):
    """Calculate the projection of the template in f, fdot
    etc. space from the apparent f, fdot etc. and the
    position and search Tspan"""
    freq = float(loudestLine[0])
    delta = float(loudestLine[2])
    f1dot = float(loudestLine[3])
    f2dot = float(loudestLine[4])
    f3dot = float(loudestLine[5])
    dfSpin = Tspan*(f3dot/6 + f2dot/2 + f1dot)
    # extra frequency band due to sidereal Doppler motion
    # constant is (earth radius)*(sidereal angular frequency)/c
    # assume detector is at equator (cos latitude = 1)
    dfDoppl = 1.55e-6 * cos( delta ) * max(freq, freq + dfSpin)
    f0 = min( freq, freq + dfSpin ) - dfDoppl
    f1 = max( freq, freq + dfSpin ) + dfDoppl
    return str(f0),str(f1-f0)
############################################################





############################################################
###             End of CasACommon.py                     ###
############################################################
