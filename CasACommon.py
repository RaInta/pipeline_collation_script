#!/usr/bin/python
#
###########################################
#
# File: CasACommon.py
# Author: Ra Inta
# Description:
# Created: October 20, 2016
# Last Modified: 20180326, R.I.
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
import re
from math import ceil, cos, sin, pi
import glob
from distutils.spawn import find_executable
from subprocess import check_output
import getpass


INITIAL_DIR = os.getcwd()
SCRIPTS = os.path.dirname(sys.argv[0])
RESULTS_DIR = os.path.abspath(INITIAL_DIR)


# Attempt to infer LIGO username and account codes, otherwise prompt for them
# in this case, a file is created and read from subsequently
ATLAS_USER = getpass.getuser()
CRED_FILENAME = os.path.join(SCRIPTS, "ligo_creds.txt")

if ATLAS_USER == "ra1":
    LIGO_USER = "ra.inta"
    ACCT_TAG = "ligo.prod.o1.cw.directedisolatedother.coherent"
elif ATLAS_USER == "owen":
    LIGO_USER = "ben.owen"
    ACCT_TAG = "ligo.prod.o1.cw.directedisolated.coherent"
elif ATLAS_USER == "sano":
    LIGO_USER = "santiago.caride"
    ACCT_TAG = "ligo.prod.o1.cw.directedisolatedother.coherent"
elif ATLAS_USER == "binod.rajbhandari":
    LIGO_USER = "binod.rajbhandari"
    ACCT_TAG = "ligo.prod.o1.cw.directedisolatedother.coherent"
elif os.path.isfile(CRED_FILENAME):
    # read local file for LIGO credentials
    with open(CRED_FILENAME) as credsFile:
        userStr = re.compile('(?<=USER:)\s*\w*\.+\w+')
        tagStr = re.compile('(?<=TAG:)\s*\w*\.+\w+')
        for subStr in credsFile.readlines():
            if userStr.search(subStr):
                LIGO_USER = userStr.search(subStr).group(0).split()[0]
            if tagStr.search(subStr):
                ACCT_TAG = tagStr.search( subStr).group(0).split()[0]
else:
    print("Your LIGO username has not been inferred.")
    LIGO_USER = raw_input("Please enter your albert.einstein username: ")
    ACCT_TAG = raw_input("Please enter your group accounting tag (e.g. ligo.prod.s5.cw.directedisolated.coherent): ")
    with open(CRED_FILENAME, 'w') as credsFile:
        credsFile.write( "LIGO_USER: " + LIGO_USER + "\n")
        credsFile.write( "ACCT_TAG: " + ACCT_TAG )




##################################################


##################################################
## Some more search-wide constants       #########
##################################################

FSTAT_DTERMS = 8  # Number of terms in Dirichlet kernel (cut from default 16 to match with SSE2 fstat code)
UL_HIST_BINWIDTH = 10

# Files and paths
#
SRCH_HISTOGRAM = os.path.join(RESULTS_DIR, 'search_histogram.txt')
UL_MISM_HISTOGRAM = os.path.join(SCRIPTS, 'upper_limit_mismatch_histogram.txt')

SETUP_XML = os.path.join(RESULTS_DIR, 'search_setup.xml')
SFT_STRETCH_XML = os.path.join(RESULTS_DIR, 'optimal_sft_stretch.xml')
COMP_COST_XML = os.path.join(RESULTS_DIR, 'compute_cost.xml')
SFT_DB = os.path.join(RESULTS_DIR, 'sft_database.xml.bz2')
DENSITY_DB = os.path.join(RESULTS_DIR, 'template_density.xml')
VETO_BANDS_DB = os.path.join(RESULTS_DIR, 'veto_bands.xml')
SRCH_BANDS_DB = os.path.join(RESULTS_DIR, 'search_bands.xml.bz2')
UL_BANDS_DB = os.path.join(RESULTS_DIR, 'upper_limit_bands.xml')

# Directory functions
def local_dir(path):
    """docstring for local_dir"""
    return os.path.join(RESULTS_DIR, 'local', path)

def global_dir(path):
    """docstring for global_dir"""
    return os.path.join(RESULTS_DIR, 'global', path)

def lal_data_path(path=''):
    """docstring for lal_data_path"""
    return os.path.join(SCRIPTS, "production", "share", "lalpulsar", path)

# Make executable path a list
EXEC_PATH = [ os.path.abspath(x) for x in [SCRIPTS, 'production/bin'] ]

##################################################


# Print hostname if running as Condor job
for envKey in os.environ:
    if re.match(".*_CONDOR", envKey):
        print("%%%%%%%%%%\n Running under Condor on " \
              + check_output("hostname --short", shell=True) \
              + "%%%%%%%%%%\n")
        break



############################################################
###   Check no files match the fileName  ###################
############################################################
def assert_no_files_matching( fileNameList ):
    """Checks that there are no files matching the (globbed list) input"""
    #TODO refactor this as try... except
    for fileGlobIdx in fileNameList:
        for fileIdx in glob.glob( fileGlobIdx ):
            if os.path.isfile( fileIdx ):
                print(os.getcwd() + " is not empty of files matching " + fileGlobIdx )
                break
############################################################

# Assert some things are strictly positive
def assert_strictly_positive(num, numIdx):
    """docstring for assert_strictly_positive"""
    assert num > 0, numIdx + " must be strictly positive"

# Assert some things are positive
def assert_positive(num, numIdx):
    """docstring for assert_positive"""
    assert num >= 0, numIdx + " must be positive"

# Assert some things are strictly less than other things
def assert_strictly_less(num1, num2):
    """docstring for assert_strictly_less"""
    assert num1 < num2, str(num1) + " must be strictly less than " + str(num2)

# Assert some things are less than other things
def assert_less():
    """docstring for assert_less"""
    assert num1 <= num2, str(num1) + " must be less than " + str(num2)

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
# Open possibly BZ2 zipped files
############################################################

def openFile(fileName):
    if fileName[-4:] == '.bz2':
        return bz2.BZ2File(fileName, "r")
    else:
        return open(fileName, "r")

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
        # Default threading, RAM and disk requests
        self.request_cpus = '1'
        self.request_memory = '8 GB'
        self.request_disk = '16 GB'
        # TODO Let's get rid of these globals --- maybe put in search_setup.xml
        global LIGO_USER
        global ACCT_TAG
        global INITIAL_DIR
        global SCRIPTS

    def write(self):
        """Writes out the Condor submit file"""
        # The boilerplate stuff
        header_str = "# global options set by " + __file__ + "\n\n"
        header_str += "universe = vanilla\n"
        header_str += "initialdir = " + os.getcwd() + "\n"
        header_str += "getenv = true\n"
        header_str += "notification = never\n"
        header_str += "log = condor.log\n"
        header_str += "request_cpus = " + self.request_cpus + "\n"
        header_str += "request_memory = " + self.request_memory + "\n"
        header_str += "request_disk = " + self.request_disk + "\n"
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
        with open(self.fileName, 'w') as subFileName:
            subFileName.write( header_str + body_str )

    def write_dict(self):
        """Writes out the Condor submit job as a (htcondor compliant) dictionary"""
        # The boilerplate stuff
        #header_str = "# global options set by " + __file__ + "\n\n"
        sub_dict = {
        "universe" : "vanilla",
        "initialdir" : INITIAL_DIR,
        "getenv" : "true",
        "notification" : "never",
        "log" : "condor.log",
        "executable": os.path.join(SCRIPTS, self.executable),
        "accounting_group" : ACCT_TAG,
        "accounting_group_user" : LIGO_USER,
        "output" : self.output,
        "error" : self.error,
        "request_cpus" : self.request_cpus,
        "request_memory" : self.request_memory,
        "request_disk" : self.request_disk,
        "queue " : str(self.queue)
        }
        if len(self.args):
            arg_str = ""
            arg_str += " ".join([str(x) for x in self.args] )
            sub_dict["arguments"] = arg_str
        # Job specific options
        #body_str = "# this job's options and commands\n\n"
        # TODO can you add the queue as a htcondor.Submit() object property?
        #body_str += "queue " + str(self.queue) + "\n"
        return sub_dict


############################################################



############################################################
###    Read setup XMLs into an object       ################
############################################################

class SetupXML:
    def __init__(self, fileName='search_setup.xml'):
        """docstring for SetupXML

        For example:
            searchSetup = SetupXML('search_setup.xml')
            searchSetup.create("setup")
            searchSetup.add("search", "band", "200")
            searchSetup.add("upper_limit", "band", "1.0")
            searchSetup.write()

        will write out the search_setup.xml, populated with the appropriate
        values:
            searchSetup.read("setup", "upper_limit", "band")
            >> "1.0"

            """
        # Get parameters from search_setup.xml
        self.fileName = fileName
        self.setup_tree = None
    def read(self, mainTag, innerTag, contentsTag):
        """docstring for read
        Reads in first and second level tags (innerTag and contentsTag), possibly
        as a list. Returns the value of the element.
        For example:
            setupObj = SetupXML('search_setup.xml')
            [distance, moment_of_inertia] = setupObj.read(['target','target'],['distance','moment_of_inertia'])
            print("Dist: " + str( distance ) + "\nI_zz: " + str( moment_of_inertia) )
            Dist: 1.0799880875e+20
            I_zz: 1e+38
            """
        self.mainTag = mainTag
        self.innerTag = innerTag
        self.contentsTag = contentsTag
        tag_contents = []
        # This will open bz2 files too, as long as fileName has a .bz2 extension
        with openFile( self.fileName ) as xmlContents:
            xml_root = ET.parse( xmlContents ).getroot()
        for main_tag in xml_root.iter( self.mainTag):
            for search_params in main_tag.iter( self.innerTag ):
                for contentsIdx in search_params.iter( self.contentsTag ):
                    #tag_contents.append( float( search_params.find( contentsIdx ).text ) )
                    tag_contents.append(  contentsIdx.text )
                    #tag_contents = contentsIdx.text
        return tag_contents
    def create(self, mainTag, comment=""):
        # TODO doesn't parse comment properly yet.
        """
        Create an XML in a format consumable by this pipeline.
        e.g.
            setupObj = SetupXML("search_setup.xml")
            setupObj.create("setup")
        """
        self.mainTag = mainTag
        self.comment = comment
        self.setup_tree = ET.Element(self.mainTag)
        self.setup_tree.Comment = self.comment
        #ET.Comment = self.comment
    def add(self, innerTag, contentsTag, valueText):
        """
        Append values to a setup XML already created using .create() method
        For example:
            setupObj = SetupXML("search_setup.xml")
            setupObj.create("setup")
            setupObj.add("target","distance", "1.0799880875e+20")
            """
        self.innerTag = innerTag
        self.contentsTag = contentsTag
        self.valueText = valueText
        alreadyExists = self.setup_tree.find( self.innerTag )
        if alreadyExists:
            inside = ET.SubElement( alreadyExists, self.contentsTag)
            inside.text = self.valueText
        else:
            inner = ET.SubElement(self.setup_tree, self.innerTag)
            inside = ET.SubElement(inner, self.contentsTag)
            inside.text = self.valueText
    def write(self):
        # TODO finish this, as it adds extra tags because it just reads in lists
        # for tags.
        """docstring for write
        Generates an ElementTree populated by the list elements innerTag and contentsTag"""
        #indent(setup_root)
        #setup_xml = ET.ElementTree(setup_root)
        #setup_xml.write(self.fileName + '.xml', xml_declaration=True, encoding='UTF-8', method='html')
        indent(self.setup_tree)
        outside_xml = ET.ElementTree( self.setup_tree )
        outside_xml.write(self.fileName, xml_declaration=True, encoding='UTF-8', method='xml')


############################################################

# Check that directory exists and is accessible
def check_directory(dir):
    """docstring for check_directory"""
    if not os.path.exists(dir):
        # Note: at least the version of OS I'm using, symlinks are also returned
        # as True with os.path.exists()
        ## It might be a symlink!
        #if os.path.islink(dir):
        #    dir = os.readlink(dir)
        os.system("mkdir -p " + dir + " >/dev/null 2>&1")
    # return Boolean if we have read, write and execute access to the directory
    return os.access(dir, os.F_OK) & os.access(dir, os.X_OK) & os.access(dir, os.R_OK) & os.access(dir, os.W_OK)

# Assert directory exists
def assert_directory(dir):
    """Asserts if dir is a directory created by check_directory()"""
    assert check_directory(dir), dir + "does not exist or is inaccessible"

############################################################
##  Check and validate SFTs    #############################
############################################################

# Check for the local, and then global, SFT directory
# TODO test this...
def get_check_local_sfts(jobNum):
    """docstring for get_check_local_sfts"""
    pass
    old_cwd = os.getcwd()
    localdir = local_dir('sfts')
    globaldir = global_dir('sfts')
    assert_directory(globaldir)
    print("----------\n")
    # Try to use local SFTs
    if check_directory(local):
        # Master checksum file
        global_checksum = os.path.join(globaldir, 'md5sums.txt')
        # local checksum files
        local_job_checksum = os.path.join(localdir, "md5sums.txt." + str(jobNum))
        local_checksums = os.path.join(localdir, "md5sums.txt.*")
        # Check if local checksum matches master checksum
        def check(local_checksum):
            for local_check in local_checksum:
                if os.path.isfile(check):
                    if check_output("diff --text --brief" + local_check + " " + global_checksum + " >/dev/null"):
                        print("SFT checksum " + local_check + " matches " + global_checksum + "\n")
                        sfts = local
                        print("Using local SFTs: " + sfts + "\n")
                        break
                    else:
                        print("SFT checksum " + local_check + " does not match " + global_checksum + "\n")
                        sys.exit(0)
                else:
                    print("SFT checksum " + local_check + " does not exist\n")
                    sys.exit(0)
        # Try all existing local checksum files
        for local_checksum in glob.glob(local_checksums):
            if check(local_checksum):
                break

        # If no valid checksums files were found
        if not 'sfts' in locals():
            # Create a checksum file
            print("Computing SFT checksum " + local_job_checksum + " ...")
            try:
                os.system("rm -f " + local_job_checksum)
            except OSError, e:
                print("Couldn't remove " + local_job_checksum + "!")
                raise e
            os.chdir(localdir)
            for sft in sorted(glob.glob('*.sft'), key=str.lower):
                try:
                    os.system("md5sum --binary " + sft + " >> " + local_job_checksum)
                except OSError, e:
                    print("md5sum failed on " + sft + "!")
                    raise e
            os.chdir(old_cwd)
            print(" done\n")
            # Try this checksum file
            check(local_job_checksum)
    # Only the global SFTs are available, use as a last resort
    if not 'sfts' in locals():
        sfts = globaldir
        print("Using global SFTs: " + sfts + "\n")

    print("----------\n")
    return sfts



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
    dfSpin = (((f3dot/6)*Tspan + f2dot/2)*Tspan + f1dot)*Tspan
    # extra frequency band due to sidereal Doppler motion
    # constant is (earth radius)*(sidereal angular frequency)/c
    # assume detector is at equator (cos latitude = 1)
    dfDoppl = 1.55e-6 * cos( delta ) * max(freq, freq + dfSpin)
    f0 = min( freq, freq + dfSpin ) - dfDoppl
    f1 = max( freq, freq + dfSpin ) + dfDoppl
    return f0, f1 - f0

############################################################


def run_exec(cmd, verbose=True):
    """Execute OS functions in a consistent way.
    This is so e.g. flags can be passed easily from each
    search step. It requires the subprocess module, which
    takes successive flags as a list."""
    cmdStr = ""
    for cmds in cmd:
        cmdStr += cmds
    if verbose:
        printStr = ""
        printStr += "##########\n"
        printStr += "Path: " + os.path.abspath(cmd[0]) + "\n"
        printStr += "Executing" + cmdStr + "\n"
        printStr += "##########\n"
        printStr += "##########\n"
    # Test if the command is executable
    if not find_executable(cmd[0]):
        print("Command " + cmd[0] + " is not a valid executable.")
        exit(0)
    # Try to run it
    try:
        subprocess.call(cmd)
        print(printStr)
    except OSError, e:
        raise cmdStr + " failed to run ", e


############################################################
# Calculate orbital Doppler factor
############################################################
# NOTE most of this is archived currently; the first section is
# calculating the Doppler shift a priori and has been replaced... for now.

## Astrophysical constants
#W_YEAR = 2*pi/(365.25*24*3600)  # Angular velocity of the Earth around the Sun
#W_DAY = 2*pi/(23.9344699*3600)  # Angular velocity of the Earth's own rotation (sidereal)
#R_AU = 149597870700  # 1 Au (m)
#R_EARTH = 6371000  # Earth's mean radius (m)
## Calculate Earth's tangential orbital velocity:
#c = 299792458.0  # m/s
#V_YEAR = W_YEAR*R_AU/c  # Dimensionless beta of Earth's solar orbital speed
#V_DAY = W_DAY*R_EARTH/c  # Dimensionless beta of Earth's own rotational speed
#
## TODO use ephemereides files like sane people do (mind you, the current method
## saves reading a file each time)
## For now, let's hard-code the 2017 Vernal Equinox:
## Monday, March 20, 2017, 10:29 UTC
#T_EQUINOX = 1174040956
#
#def deg2rad(degrees, minutes=0.0, seconds=0.0):
#    """Converts from degrees, minutes, seconds to radians"""
#    return (pi/180)*(degrees+(minutes/60.0)+(seconds/3600.0))
#
#def epsilon(t):
#    """Calculates the Earth's orbital inclination to the ecliptic at
#    time t (GPS).
#
#    This uses JPL's updated ephemerides:
#    epsilon = 23deg 26m 21.406s -46.836769s * T -0.0001831s * pow(T, 2)
#    +0.00200340s * pow(T, 3) -0.576e-6s * pow(T, 4) -4.34e-8s * pow(T, 5)
#    Where T is the time, in Julian centuries, since the J2000 epoch.
#    [Astronomical Almanac 2010, p. B52]
#
#    GPS time for J2000 epoch:
#    January 1, 2000, 11:58:55.816 UTC"""
#    t_J2000 = 630763148
#    T = (t - t_J2000)/(100*365.25*24*3600)
#    return deg2rad(23, 26, 21.406) - deg2rad(0, 0, 46.836769)*T - deg2rad(0, 0, 0.0001831)*pow(T, 2) + deg2rad(0, 0, 0.00200340)*pow(T, 3) \
#           - deg2rad(0, 0, 0.576e-6)*pow(T, 4) - deg2rad(0, 0, 4.34e-8) * pow(T, 5)
#
#def orbitalDoppler(alpha, delta, t):
#    """A function to calculate the orbital component of the Doppler shift due to Earth's
#    orbital motion around the Sun. It requires Right Ascension (alpha), declination (delta),
#    and GPS time (t). It also requires the time of the (Northern Hemisphere) Vernal Equinox,
#    and the inclination to the ecliptic (epsilon) although these are referenced internally.
#    The Doppler factor is returned."""
#    epsilon_t = epsilon(t)
#    return 1 -V_YEAR*(cos(delta)*sin(alpha)*cos(epsilon_t) + sin(delta)*sin(epsilon_t))*cos(W_YEAR*(t - T_EQUINOX)) \
#           +V_YEAR*cos(delta)*cos(alpha)*cos(W_YEAR*(t - T_EQUINOX) )
#
#def siderealDoppler(alpha, delta, t):
#    """Calculates Doppler shift factor as a result of Earth's sidereal motion"""
#    # TODO Don't use this yet---needs fixing
#    epsilon_t = epsilon(t)
#    # \del(v.n)
#    #return -V_DAY*W_DAY*delta_t *(cos(delta)*sin(alpha)*cos(epsilon_t) + sin(delta)*sin(epsilon_t))*cos(W_DAY*(t-T_EQUINOX)) \
#    #        -V_DAY*W_DAY*delta_t*cos(delta)*cos(alpha)*sin(W_DAY*(t-T_EQUINOX))
#    return 1 -V_DAY*W_DAY*(cos(delta)*sin(alpha)*cos(epsilon_t) + sin(delta)*sin(epsilon_t))*cos(W_DAY*(t-T_EQUINOX)) \
#            -V_DAY*W_DAY*cos(delta)*cos(alpha)*sin(W_DAY*(t-T_EQUINOX))

def orbitalDoppler(alpha, delta, t):
    """A function to calculate the orbital component of the Doppler shift due to Earth's
    orbital motion around the Sun. It requires Right Ascension (alpha), declination (delta),
    and GPS time (t). It also requires the time of the (Northern Hemisphere) Vernal Equinox,
    and the inclination to the ecliptic (epsilon) although these are referenced internally.
    The Doppler factor is returned."""
    # Call the lalapps function
    execPath = os.path.abspath(os.path.join(SCRIPTS, "production", "bin", "lalapps_PrintDetectorState"))
    ephemPath = os.path.join(SCRIPTS, "production", "share", "lalpulsar")
    DF_H1 = check_output([execPath, "-I", "H1", "-a", str(alpha),"-d", str(delta), "-t", str(t), "-y", "00-19-DE405", "-E", ephemPath])
    srchStr = re.compile(r'dtSSB/dtDet - 1 = [-+]?[0-9]+.[0-9]+e[-+]?[0-9]+')
    srchNum = re.compile(r'[-+]?[0-9]+.[0-9]+e[-+]?[0-9]+')
    rateStr = srchNum.search(srchStr.search(DF_H1).group()).group()
    return 1 + float(rateStr)

############################################################
###             End of CasACommon.py                     ###
############################################################
