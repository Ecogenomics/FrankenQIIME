#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jens Reeder","Dan Knights", "Antonio Gonzalez Pena",
    "Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

from qiime.util import make_option
from os import access, X_OK, R_OK, W_OK, getenv, environ
from os.path import isdir, exists, split
from sys import platform, version as python_version, executable
from shutil import rmtree
from subprocess import Popen, PIPE, STDOUT

from cogent.util.unit_test import TestCase, main as test_main
from cogent.util.misc import app_path, get_random_directory_name
from cogent.app.util import ApplicationNotFoundError, ApplicationError

from qiime.parse import parse_qiime_config_file
from qiime.util import (load_qiime_config, 
                        get_qiime_project_dir, 
                        parse_command_line_parameters,
                        get_qiime_library_version,
                        get_rdp_jarpath)
from qiime.denoiser.utils import check_flowgram_ali_exe

from cogent import __version__ as pycogent_lib_version
from numpy import __version__ as numpy_lib_version
try:
    from matplotlib import __version__ as matplotlib_lib_version
except ImportError:
    matplotlib_lib_version = "Not installed."
try:
    from pynast import __version__ as pynast_lib_version
except ImportError:
    pynast_lib_version = "Not installed."

script_info = {}
script_info['brief_description']= """Print out the qiime config settings."""
script_info['script_description'] = """A simple scripts that prints out the qiime config settings and does some sanity checks."""
script_info['script_usage']=[]
script_info['script_usage'].append(
    ("Example 1","""Print qiime config settings:""","""print_qiime_config.py"""))
script_info['script_usage'].append(
    ("Example 2","""Print and check qiime config settings for sanity:""",
     """print_qiime_config.py -t"""))

script_info['output_description'] = """This prints the qiime_config to stdout."""
script_info['version'] = __version__
script_info['help_on_no_arguments'] = False
script_info['required_options']=[]
script_info['optional_options'] = [\
    make_option('-t','--test', action='store_true',
                dest='test', default = False,
                help='Test the qiime config for sanity '
                +'[default: %default]')]

class Qiime_config(TestCase):
    
    def setUp(self):
        self.config = load_qiime_config()
    
    def test_python_exe_fp(self):
        """python_exe_fp is set to a working python env"""
        
        python = self.config["python_exe_fp"]
        command = "%s --version" % python
        proc = Popen(command,shell=True,universal_newlines=True,\
                         stdout=PIPE,stderr=STDOUT)
        #Check if callable
        if proc.wait() !=0:
            self.fail("Calling python failed. Check you python_exe_fp:%s" %python)
        
        # Does it give its version string?
        out_string = proc.stdout.read()
        if not out_string:
            self.fail("Something is wrong with your python\n." \
                          +" Check you python_exe_fp:%s" %python)
    
    def test_cluster_jobs_fp(self):
        """cluster_jobs_fp is set to a valid path and is executable"""       
        
        test_qiime_config_variable("cluster_jobs_fp", self.config, self, X_OK,)
   
    def test_blastmat_dir(self):
        """blastmat_dir is set to a valid path."""
        
        test_qiime_config_variable("blastmat_dir", self.config, self)
        
    def test_blastall_fp(self):
        """blastall_fp is set to a valid path"""
        
        blastall = self.config["blastall_fp"]
        if not self.config["blastall_fp"].startswith("/"):
            #path is relative, figure out absolute path
            blast_all = app_path(blastall)
            if not blast_all:
                raise ApplicationNotFoundError("blastall_fp set to %s, but is not in your PATH. Either use an absolute path to or put it in your PATH." % blastall)
            self.config["blastall_fp"] = blast_all

        test_qiime_config_variable("blastall_fp", self.config, self, X_OK)

    def test_ampliconnoise_install(self):
        """ AmpliconNoise install looks sane."""
        url = "http://www.qiime.org/install/install.html#ampliconnoise-install"
        
        pyro_lookup_file = getenv('PYRO_LOOKUP_FILE')
        self.assertTrue(pyro_lookup_file != None,
         "$PYRO_LOOKUP_FILE variable is not set. See %s for help." % url)
        self.assertTrue(exists(pyro_lookup_file),
         "$PYRO_LOOKUP_FILE variable is not set to an existing filepath.")
         
        seq_lookup_file = getenv('SEQ_LOOKUP_FILE')
        self.assertTrue(seq_lookup_file != None,
         "$SEQ_LOOKUP_FILE variable is not set. See %s for help." % url)
        self.assertTrue(exists(seq_lookup_file),
         "$SEQ_LOOKUP_FILE variable is not set to an existing filepath.")
         
        self.assertTrue(app_path("SplitKeys.pl"),
         "AmpliconNoise Scripts directory doesn't seem to be in $PATH."+\
         " See %s for help." % url)
         
        self.assertTrue(app_path("Perseus"),
         "AmpliconNoise bin directory doesn't seem to be in $PATH."+\
         " See %s for help." % url)
        
    def test_pynast_template_alignment_fp(self):
        """pynast_template_alignment, if set, is set to a valid path"""
            
        test_qiime_config_variable("pynast_template_alignment_fp",
                                   self.config, self)
            
    def test_pynast_template_alignment_blastdb_fp(self):
        """pynast_template_alignment_blastdb, if set, is set to a valid path"""
            
        test_qiime_config_variable("pynast_template_alignment_blastdb_fp",
                                   self.config, self)
    def test_pynast_template_alignment_blastdb_fp(self):
        """pynast_template_alignment_blastdb, if set, is set to a valid path"""
        
        test_qiime_config_variable("pynast_template_alignment_blastdb_fp",
                                   self.config, self)
        
    def test_template_alignment_lanemask_fp(self):
        """template_alignment_lanemask, if set, is set to a valid path"""
            
        test_qiime_config_variable("template_alignment_lanemask_fp",
                                   self.config, self)
    
    def test_qiime_scripts_dir(self):
        """qiime_scripts_dir, if set, is set to a valid path"""

        scripts_dir = self.config["qiime_scripts_dir"]
        
        if scripts_dir:
            self.assertTrue(exists(scripts_dir),
                            "qiime_scripts_dir does not exist: %s" % scripts_dir)
            self.assertTrue(isdir(scripts_dir),
                            "qiime_scripts_dir is not a directory: %s" % scripts_dir)
        else:
            pass
            #self.fail("scripts_dir is not set.")

    def test_temp_dir(self):
        """temp_dir, if set, is set to a valid path"""

        temp_dir = self.config["temp_dir"]
        
        if temp_dir:
            self.assertTrue(exists(temp_dir),
                            "temp_dir does not exist: %s" % temp_dir)
            self.assertTrue(isdir(temp_dir),
                            "temp_dir is not a directory: %s" % temp_dir)
        else:
            pass
            #self.fail("temp_dir is not set.")

    def test_working_dir(self):
        """working_dir, if set, is set to a valid path"""

        working_dir = self.config["working_dir"]
        
        if working_dir:
            self.assertTrue(exists(working_dir), 
                            "working dir does not exist: %s" % working_dir)
            self.assertTrue(isdir(working_dir),
                            "working_dir is not a directory: %s" % working_dir)        
            self.assertTrue(access(working_dir, W_OK),
                            "working_dir not writable: %s" % working_dir)
        else:
            pass
            #self.fail("working_dir is not set.")
        
    #we are not testing these values from the qiime_config:
    # jobs_to_start   1
    # seconds_to_sleep        60


    def test_for_obsolete_values(self):
        """local qiime_config has no extra params"""
        
        qiime_project_dir = get_qiime_project_dir()
        orig_config = parse_qiime_config_file(open(qiime_project_dir +
                                             '/qiime/support_files/qiime_config'))
        
        #check the env qiime_config
        qiime_config_env_filepath = getenv('QIIME_CONFIG_FP')
        if qiime_config_env_filepath:
            qiime_config_via_env= parse_qiime_config_file(open(qiime_config_env_filepath))
            extra_vals = []
            for key in qiime_config_via_env:
                if key not in orig_config:
                    extra_vals.append(key)
            if extra_vals:
                self.fail("The qiime_config file set via QIIME_CONFIG_FP"+
                          "enviroment variable contains obsolete parameters:\n"+
                          ", ".join(extra_vals))
        # check the qiime_config in $HOME/.qiime_config
        home_dir = getenv('HOME')        
        if (exists(home_dir+"/.qiime_config")):
            qiime_config_home = parse_qiime_config_file(open(home_dir+"/.qiime_config"))
            extra_vals = []
            for key in qiime_config_home:
                if key not in orig_config:
                    extra_vals.append(key)
            if extra_vals:
                self.fail("The .qiime_config in your HOME contains obsolete "+
                          "parameters:\n" + ", ".join(extra_vals))

    def test_chimeraSlayer_install(self):
        """no obvious problems with ChimeraSlayer install """

        #The ChimerSalyer app requires that all its components are installed
        # relative to the main program ChimeraSlayer.pl.
        # We therefore check that at least one the files is there.
        # However, if the directory structure of ChimeraSlayer changes, this test will most
        # likely fail as well and need to be updated.
        # Tested with the version of microbiomeutil_2010-04-29

        chim_slay = app_path("ChimeraSlayer.pl")
        self.assertTrue(chim_slay,"ChimeraSlayer was not found in your $PATH")
        dir, app_name = split(chim_slay)
        self.assertTrue(exists(dir+"/ChimeraParentSelector/chimeraParentSelector.pl"),
         "ChimeraSlayer depends on external files in directoryies relative to its "
         "install directory. Thesedo not appear to be present.")

    def test_uclust_supported_version(self):
        """uclust is in path and version is supported """
        acceptable_version = (1,2,22)
        self.assertTrue(app_path('uclust'),
         "uclust not found. This may or may not be a problem depending on "+\
         "which components of QIIME you plan to use.")
        command = 'uclust --version'
        proc = Popen(command,shell=True,universal_newlines=True,\
                         stdout=PIPE,stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip().split('v')[-1].strip('q')
        try:
            version = tuple(map(int,version_string.split('.')))
            pass_test = version == acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,\
         "Unsupported uclust version. %s is required, but running %s." \
         % ('.'.join(map(str,acceptable_version)), version_string))

    def test_python_supported_version(self):
        """python is in path and version is supported """
        min_acceptable_version = (2,6,0)
        min_unacceptable_version = (2,7,0)
        command = 'python --version'
        proc = Popen(command,shell=True,universal_newlines=True,\
                         stdout=PIPE,stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip().split('Python')[-1].strip()
        try:
            version = tuple(map(int,version_string.split('.')))
            if len(version) == 2:
                version = tuple(version[0],version[1],0)
            pass_test = (version >= min_acceptable_version and version < min_unacceptable_version)
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,\
         "Unsupported python version. Must be >= %s and < %s , but running %s." \
         % ('.'.join(map(str,min_acceptable_version)),
            '.'.join(map(str,min_unacceptable_version)),
            version_string))

    def test_numpy_suported_version(self):
        """numpy version is supported """
        min_acceptable_version = (1,3,0)
        min_unacceptable_version = (1,5,1)
        try:
            from numpy import __version__ as numpy_lib_version
            version = tuple(map(int,numpy_lib_version.split('.')))
            pass_test = (version >= min_acceptable_version and version < min_unacceptable_version)
            version_string = str(numpy_lib_version)
        except ImportError:
            pass_test = False
            version_string = "Not installed"
        self.assertTrue(pass_test,\
         "Unsupported numpy version. Must be >= %s and < %s , but running %s." \
         % ('.'.join(map(str,min_acceptable_version)),
            '.'.join(map(str,min_unacceptable_version)),
            version_string))

    def test_matplotlib_suported_version(self):
        """maptplotlib version is supported """
        min_acceptable_version = (0,98,5,3)
        min_unacceptable_version = (0,98,5,4)
        try:
            from matplotlib import __version__ as matplotlib_lib_version
            version = tuple(map(int,matplotlib_lib_version.split('.')))
            pass_test = (version >= min_acceptable_version and version < min_unacceptable_version)
            version_string = str(matplotlib_lib_version)
        except ImportError:
            pass_test = False
            version_string = "Not installed"
        self.assertTrue(pass_test,\
         "Unsupported matplotlib version. Must be >= %s and < %s , but running %s." \
         % ('.'.join(map(str,min_acceptable_version)),
            '.'.join(map(str,min_unacceptable_version)),
            version_string))
            
    def test_pynast_suported_version(self):
        """pynast version is supported """
        min_acceptable_version = (1,1)
        min_unacceptable_version = (1,2)
        try:
            from pynast import __version__ as pynast_lib_version
            version = pynast_lib_version.split('.')
            if version[-1][-4:]=='-dev':
                 version[-1] = version[-1][:-4]
            version = tuple(map(int,version))
            pass_test = (version >= min_acceptable_version and version < min_unacceptable_version)
            version_string = str(pynast_lib_version)
        except ImportError:
            pass_test = False
            version_string = "Not installed"
        self.assertTrue(pass_test,\
         "Unsupported matplotlib version. Must be >= %s and < %s , but running %s." \
         % ('.'.join(map(str,min_acceptable_version)),
            '.'.join(map(str,min_unacceptable_version)),
            version_string))

    def test_blast_supported_version(self):
        """blast is in path and version is supported """
        acceptable_version = (2,2,22)
        self.assertTrue(app_path('blastall'),
         "blast not found. This may or may not be a problem depending on "+\
         "which components of QIIME you plan to use.")
        command = 'blastall | grep blastall'
        proc = Popen(command,shell=True,universal_newlines=True,\
                         stdout=PIPE,stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip().split(' ')[1].strip()
        try:
            version = tuple(map(int,version_string.split('.')))
            pass_test = version == acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,\
         "Unsupported blast version. %s is required, but running %s." \
         % ('.'.join(map(str,acceptable_version)), version_string))
         
    def test_FastTree_supported_version(self):
        """FastTree is in path and version is supported """
        acceptable_version = (2,1,3)
        self.assertTrue(app_path('FastTree'),
         "FastTree not found. This may or may not be a problem depending on "+\
         "which components of QIIME you plan to use.")
        command = "FastTree 2>&1 > /dev/null | grep version"
        proc = Popen(command,shell=True,universal_newlines=True,\
                         stdout=PIPE,stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip().split(' ')[4].strip()
        try:
            version = tuple(map(int,version_string.split('.')))
            pass_test = version == acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,\
         "Unsupported FastTree version. %s is required, but running %s." \
         % ('.'.join(map(str,acceptable_version)), version_string))
    
    def test_cdbtools_supported_version(self):
        """cdbtools is in path and version is supported """
        acceptable_version = (0,99)
        self.assertTrue(app_path('cdbfasta'),
         "cdbtools not found. This may or may not be a problem depending on "+\
         "which components of QIIME you plan to use.")
        command = "cdbfasta -v"
        proc = Popen(command,shell=True,universal_newlines=True,\
                         stdout=PIPE,stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip().split(' ')[2].strip()
        try:
            version = tuple(map(int,version_string.split('.')))
            pass_test = version == acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,\
         "Unsupported cdbtools version. %s is required, but running %s." \
         % ('.'.join(map(str,acceptable_version)), version_string))
        
    def test_INFERNAL_supported_version(self):
        """INFERNAL is in path and version is supported """
        acceptable_version = (1,0,2)
        self.assertTrue(app_path('cmbuild'),
         "Infernal not found. This may or may not be a problem depending on "+\
         "which components of QIIME you plan to use.")
        command = "cmbuild -h | grep INF"
        proc = Popen(command,shell=True,universal_newlines=True,\
                         stdout=PIPE,stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip().split(' ')[2].strip()
        try:
            version = tuple(map(int,version_string.split('.')))
            pass_test = version == acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,\
         "Unsupported INFERNAL version. %s is required, but running %s." \
         % ('.'.join(map(str,acceptable_version)), version_string))
        
    def test_muscle_supported_version(self):
        """muscle is in path and version is supported """
        acceptable_version = (3,6)
        self.assertTrue(app_path('muscle'),
         "muscle not found. This may or may not be a problem depending on "+\
         "which components of QIIME you plan to use.")
        command = "muscle -version"
        proc = Popen(command,shell=True,universal_newlines=True,\
                         stdout=PIPE,stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip().split(' ')[1].strip('v')
        try:
            version = tuple(map(int,version_string.split('.')))
            pass_test = version == acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,\
         "Unsupported muscle version. %s is required, but running %s." \
         % ('.'.join(map(str,acceptable_version)), version_string))
        
    def test_mothur_supported_version(self):
        """mothur is in path and version is supported """
        acceptable_version = (1,6,0)
        self.assertTrue(app_path('mothur'),
         "mothur not found. This may or may not be a problem depending on "+\
         "which components of QIIME you plan to use.")
        # mothur creates a log file in cwd, so create a tmp and cd there first
        tmp_dir = get_random_directory_name(output_dir='/tmp/')
        command = "cd %s ; mothur -v | grep ^mothur" % tmp_dir
        proc = Popen(command,shell=True,universal_newlines=True,\
                         stdout=PIPE,stderr=STDOUT)
        stdout, stderr = proc.communicate()
        version_string = stdout.strip().split(' ')[1].strip('v.')
        try:
            version = tuple(map(int,version_string.split('.')))
            pass_test = version == acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,\
         "Unsupported mothur version. %s is required, but running %s." \
         % ('.'.join(map(str,acceptable_version)), version_string))
        
        # remove the directory and the log file
        rmtree(tmp_dir)
         
    def test_denoiser_supported_version(self):
        """denoiser aligner is ready to use """

        pass_test = True
        try:
            check_flowgram_ali_exe()
        except (ApplicationNotFoundError, ApplicationError):
            pass_test = False
            
        self.assertTrue(pass_test, "Denoiser flowgram aligner not found or not executable."+\
                            "This may or may not be a problem depending on "+\
                            "which components of QIIME you plan to use.")
        
    def test_raxmlHPC_supported_version(self):
        """raxmlHPC is in path and version is supported """
        acceptable_version = (7,0,3)
        self.assertTrue(app_path('raxmlHPC'),
         "raxmlHPC not found. This may or may not be a problem depending on "+\
         "which components of QIIME you plan to use.")
        command = "raxmlHPC -v | grep version"
        proc = Popen(command,shell=True,universal_newlines=True,\
                         stdout=PIPE,stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip().split(' ')[4].strip()
        try:
            version = tuple(map(int,version_string.split('.')))
            pass_test = version == acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,\
         "Unsupported raxmlHPC version. %s is required, but running %s." \
         % ('.'.join(map(str,acceptable_version)), version_string))
        
    def test_clearcut_supported_version(self):
        """clearcut is in path and version is supported """
        acceptable_version = (1,0,9)
        self.assertTrue(app_path('clearcut'),
         "clearcut not found. This may or may not be a problem depending on "+\
         "which components of QIIME you plan to use.")
        command = "clearcut -V"
        proc = Popen(command,shell=True,universal_newlines=True,\
                         stdout=PIPE,stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip().split(' ')[2].strip()
        try:
            version = tuple(map(int,version_string.split('.')))
            pass_test = version == acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,\
         "Unsupported clearcut version. %s is required, but running %s." \
         % ('.'.join(map(str,acceptable_version)), version_string))
        
    def test_cdhit_supported_version(self):
        """cd-hit is in path and version is supported """
        self.assertTrue(app_path('cd-hit'),
         "cd-hit not found. This may or may not be a problem depending on "+\
         "which components of QIIME you plan to use.")
        # cd-hit does not have a version print in their program
    
    #def test_rdp_classifier_supported_version(self):
    #    """rdp_classifier is in path and version is supported """
    #    #rdp_classifier does not have a version print in their program, not in their
    #    #command line or in their java manifest, if in the future they add it to the
    #    #manifest the command is:
    #    #unzip -c self.config['rdp_classifier_fp'] META-INF/MANIFEST.MF
    #    pass
        
    #def test_ChimeraSlayer_supported_version(self):
    #    """ChimeraSlayer is in path and version is supported """
    #    #chim_slay = app_path("ChimeraSlayer.pl")
    #    #ChimeraSlayer does not have a version print in their program
    #    pass
        
        
def test_qiime_config_variable(variable, qiime_config, test,
                               access_var=R_OK, fail_on_missing=False):
    """test if a variable is set and set to a readable path."""
    
    fp = qiime_config[variable]
    
    if not fp:
        if fail_on_missing:
            test.fail("%s not set."%variable)
        else:
            # non-essential file, so do not fail
            return
    #test if file exists    
    test.assertTrue(exists(fp), "%s set to an invalid file path: %s" %\
                        (variable,fp))

    modes = {R_OK:"readable",
             W_OK:"writable",
             X_OK:"executable"}
    #test if file readable    
    test.assertTrue(access(fp, access_var),
                    "%s is not %s: %s" % (variable, modes[access_var], fp))

if __name__ == "__main__":
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    qiime_config = load_qiime_config()
    
    rdp_jarpath = get_rdp_jarpath()
    if rdp_jarpath == None:
        rdp_version = "Not installed."
    else:
        rdp_version = split(rdp_jarpath)[1]

    system_info = [
     ("Platform", platform),
     ("Python version",python_version.replace('\n', ' ')),
     ("Python executable",executable)]
    max_len =  max([len(e[0]) for e in system_info])
    print "\nSystem information"
    print  "==================" 
    for v in system_info:
        print "%*s:\t%s" % (max_len,v[0],v[1])

    version_info = [
     ("PyCogent version", pycogent_lib_version),
     ("NumPy version", numpy_lib_version),
     ("matplotlib version", matplotlib_lib_version),
     ("QIIME library version", get_qiime_library_version()),
     ("QIIME script version", __version__),
     ("PyNAST version (if installed)", pynast_lib_version),
     ("RDP Classifier version (if installed)", rdp_version)]
    max_len =  max([len(e[0]) for e in version_info])
    print "\nDependency versions"
    print  "===================" 
    for v in version_info:
        print "%*s:\t%s" % (max_len,v[0],v[1])
    
    print "\nQIIME config values"
    print  "==================="    
    max_len =  max([len(key) for key in qiime_config])
    for key,value in  qiime_config.items():
        print "%*s:\t%s"%(max_len,key,value)

    #run the Testcase.main function to do the tests
    # need to mess with the arg string, otherwise TestCase complains
    if (opts.test):
        print "\n\nrunning checks:\n"
        test_main(argv=["","-v"])
