#!/usr/bin/env python
# File created on 20 Oct 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso","Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"
 

from cogent.util.unit_test import TestCase, main
from qiime.parallel.beta_diversity import assemble_distance_matrix
from qiime.parse import parse_distmat_to_dict
from qiime.util import get_qiime_scripts_dir
import qiime.util

import tempfile
import string
import random
import os
import shutil
import subprocess

class ParallelBetaDiversityTests(TestCase):

    def setUp(self):
        """ """
        self.dm_comp1 = dm_comp1.split('\n')
        self.dm_comp2 = dm_comp2.split('\n')
        self.dm_comp3 = dm_comp3.split('\n')
        self.dm_comp4 = dm_comp4.split('\n')
        self.expected = expected.split('\n')

        self.dirs_to_remove = []

    def tearDown(self):
        """ """

        for d in self.dirs_to_remove:
            if os.path.exists(d):
                shutil.rmtree(d)
    def test_parallel_beta_diversity(self):
        """parallel should be the same dist mtx as serial"""
        qiime_config = qiime.util.load_qiime_config()
        tempdir = qiime_config['temp_dir'] or tempfile.gettempdir() 
        # tempfile may not work on cluster, if e.g. /tmp isn't mirrored via nfs
        maindir = os.path.join(tempdir,
         ''.join(random.choice(string.ascii_letters + string.digits) \
         for x in range(10)))


        os.makedirs(maindir)
        self.dirs_to_remove.append(maindir)
        otuf = os.path.join(maindir,'otuf')
        treef = os.path.join(maindir,'treef')

        otufh = open(otuf,'w')
        otufh.write(tutorial_otu_table)
        otufh.close()

        treefh = open(treef,'w')
        treefh.write(tutorial_tree)
        treefh.close()

        scripts_dir = get_qiime_scripts_dir()
        # parallel
        cmd = scripts_dir+'/parallel_beta_diversity.py -T -O 3 '+\
         '--retain_temp_files -i %s -o %s -m unifrac -t %s' %\
         (otuf, maindir+'/para1', treef)
        # -T so doesn't return yet
        proc = subprocess.Popen(cmd,shell=True, 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        betaout, betaerr = proc.communicate()
        # first paralell version
        if betaout or betaerr:
            raise RuntimeError("parallel_beta_diversity.py should should "+\
              "not generate stdout or stderr. results:" + betaout + betaerr)

        # retain temp files doesn't matter, we just delete the folder

        # now with serial bdiv
        cmd=scripts_dir+'/beta_diversity.py -i %s -o %s -m unifrac -t %s' %\
         (otuf, maindir+'/serial1', treef)
        proc = subprocess.Popen(cmd,shell=True, 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        betaout, betaerr = proc.communicate()
        if betaout or betaerr:
            raise RuntimeError("parallel_beta_diversity.py should should "+\
              "not generate stdout or stderr. results:" + betaout + betaerr)


        serialdist =\
            parse_distmat_to_dict(open(maindir+'/serial1/unifrac_otuf','U'))

        paradist =\
            parse_distmat_to_dict(open(maindir+'/para1/unifrac_otuf','U'))

        web_res = open(maindir+'/web_res','w')

        ## use unifrac_dmtx below, from fast unifrac website march 2011
        web_res.write(unifrac_dmtx)
        web_res.close()
        unifdist = parse_distmat_to_dict(open(maindir+'/web_res','U'))
        self.assertFloatEqual(serialdist, unifdist)
        self.assertFloatEqual(paradist, unifdist)

    def test_parallel_beta_diversity2(self):
        """ extra tips in tree should not affect the unifrac dmtx"""
        qiime_config = qiime.util.load_qiime_config()
        tempdir = qiime_config['temp_dir'] or tempfile.gettempdir() 
        # tempfile may not work on cluster, if e.g. /tmp isn't mirrored via nfs
        maindir = os.path.join(tempdir,
         ''.join(random.choice(string.ascii_letters + string.digits) \
         for x in range(10)))


        os.makedirs(maindir)
        self.dirs_to_remove.append(maindir)
        otuf = os.path.join(maindir,'otuf')
        treef = os.path.join(maindir,'treef')
        trim_treef = os.path.join(maindir,'trim_treef')

        otufh = open(otuf,'w')
        otufh.write(big_otu_table)
        otufh.close()

        treefh = open(treef,'w')
        treefh.write(big_tree)
        treefh.close()

        treefh = open(trim_treef,'w')
        treefh.write(big_tree_trimmed)
        treefh.close()

        scripts_dir = get_qiime_scripts_dir()


        # parallel
        cmd = scripts_dir+'/parallel_beta_diversity.py -T -O 3 '+\
         '--retain_temp_files -i %s -o %s -m unifrac -t %s' %\
         (otuf, maindir+'/para1', treef)
        # -T so doesn't return yet
        proc = subprocess.Popen(cmd,shell=True, 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        betaout, betaerr = proc.communicate()
        if betaout or betaerr:
            raise RuntimeError("parallel_beta_diversity.py should should "+\
              "not generate stdout or stderr. results:" + betaout + betaerr)

        # parallel on trimmed_tree
        cmd = scripts_dir+'/parallel_beta_diversity.py -T -O 3 '+\
         '--retain_temp_files -i %s -o %s -m unifrac -t %s' %\
         (otuf, maindir+'/para_trim', trim_treef)
        proc = subprocess.Popen(cmd,shell=True, 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        betaout, betaerr = proc.communicate()
        if betaout or betaerr:
            raise RuntimeError("parallel_beta_diversity.py should should "+\
              "not generate stdout or stderr. results:" + betaout + betaerr)


        # retain temp files doesn't matter, we just delete the folder

        # now with serial bdiv
        cmd=scripts_dir+'/beta_diversity.py -i %s -o %s -m unifrac -t %s' %\
         (otuf, maindir+'/serial1', treef)
        proc = subprocess.Popen(cmd,shell=True, 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        betaout, betaerr = proc.communicate()
        if betaout or betaerr:
            raise RuntimeError(betaout + betaerr)


        serialdist =\
            parse_distmat_to_dict(open(maindir+'/serial1/unifrac_otuf','U'))

        paradist =\
            parse_distmat_to_dict(open(maindir+'/para1/unifrac_otuf','U'))

        paradist_trim =\
            parse_distmat_to_dict(open(maindir+'/para_trim/unifrac_otuf','U'))

        expected = open(maindir+'/expected','w')

        ## use unifrac_dmtx below, from fast unifrac website march 2011
        expected.write(big_dmtx)
        expected.close()
        unifdist = parse_distmat_to_dict(open(maindir+'/expected','U'))
        sam_keys = unifdist.keys()
        for i in range(len(sam_keys)):
            key_i = sam_keys[i]
            for j in range(i):
                key_j = sam_keys[j]
                self.assertFloatEqual(serialdist[key_i][key_j], unifdist[key_i][key_j])
                self.assertFloatEqual(paradist[key_i][key_j], unifdist[key_i][key_j])
                self.assertFloatEqual(paradist_trim[key_i][key_j], unifdist[key_i][key_j])

    def test_assemble_distance_matrix(self):
        """ assemble_distance_matrix functions as expected for a non-symmetric dm
        """
        actual = parse_distmat_to_dict(assemble_distance_matrix(
         [self.dm_comp1,self.dm_comp2,self.dm_comp3,self.dm_comp4]).split('\n'))
        exp = parse_distmat_to_dict(self.expected)
        self.assertEqual(actual,exp)
        
dm_comp1 = """	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636
PC.354	0.0	0.9999	0.637204555763	0.596124964451	0.59620370837	0.736037238732	0.790494137353	0.70551354446	0.758236466255
PC.355	0.623372333874	0.0	0.613486299912	0.631198829777	0.672502144191	0.766682588739	0.732473435165	0.738434734323	0.687363838411"""

dm_comp2 = """	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636
PC.356	0.637204555763	0.613486299912	0.0	0.691902878696	0.742922266508	0.740684448418	0.78590671873	0.723519724117	0.754301399992
PC.481	0.596124964451	0.631198829777	0.691902878696	0.0	0.668731241507	0.690956885033	0.652257003318	0.652311676944	0.661223123598"""

dm_comp3 = """	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636
PC.593	0.59620370837	0.672502144191	0.742922266508	0.668731241507	0.0	0.74048244192	0.760255173866	0.761741286862	0.751777287666
PC.607	0.736037238732	0.766682588739	0.740684448418	0.690956885033	0.74048244192	0.0	0.73093223871	0.688768250911	0.729047763792"""

dm_comp4 = """	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636
PC.634	0.790494137353	0.732473435165	0.78590671873	0.652257003318	0.760255173866	0.73093223871	0.0	0.583577161959	0.596033900207
PC.635	0.70551354446	0.738434734323	0.723519724117	0.652311676944	0.761741286862	0.688768250911	0.583577161959	0.0	0.6271222318
PC.636	0.758236466255	0.687363838411	0.754301399992	0.661223123598	0.751777287666	0.729047763792	0.596033900207	0.6271222318	0.0"""

expected = """	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636
PC.354	0.0	0.9999	0.637204555763	0.596124964451	0.59620370837	0.736037238732	0.790494137353	0.70551354446	0.758236466255
PC.355	0.623372333874	0.0	0.613486299912	0.631198829777	0.672502144191	0.766682588739	0.732473435165	0.738434734323	0.687363838411
PC.356	0.637204555763	0.613486299912	0.0	0.691902878696	0.742922266508	0.740684448418	0.78590671873	0.723519724117	0.754301399992
PC.481	0.596124964451	0.631198829777	0.691902878696	0.0	0.668731241507	0.690956885033	0.652257003318	0.652311676944	0.661223123598
PC.593	0.59620370837	0.672502144191	0.742922266508	0.668731241507	0.0	0.74048244192	0.760255173866	0.761741286862	0.751777287666
PC.607	0.736037238732	0.766682588739	0.740684448418	0.690956885033	0.74048244192	0.0	0.73093223871	0.688768250911	0.729047763792
PC.634	0.790494137353	0.732473435165	0.78590671873	0.652257003318	0.760255173866	0.73093223871	0.0	0.583577161959	0.596033900207
PC.635	0.70551354446	0.738434734323	0.723519724117	0.652311676944	0.761741286862	0.688768250911	0.583577161959	0.0	0.6271222318
PC.636	0.758236466255	0.687363838411	0.754301399992	0.661223123598	0.751777287666	0.729047763792	0.596033900207	0.6271222318	0.0"""


tutorial_tree = """(((((381:0.0213,(214:0.03728,253:0.00015)0.945:0.03224)0.763:0.00483,((269:0.02693,(231:0.00509,(105:0.01425,141:0.02641)0.846:0.01405)0.428:0.00519)0.622:0.00014,404:0.00524)0.795:0.00514)0.773:0.00508,(131:0.00518,(33:0.01631,((284:0.00828,(176:0.03098,388:0.01236)0.901:0.02175)0.885:0.01273,52:0.01046)0.743:0.00498)0.924:0.01603)0.779:0.00511)0.772:0.00014,153:0.00507)0.753:0.00602,(223:0.03237,(172:0.01733,81:0.03834)0.224:0.00414)0.845:0.01076,(136:0.00627,((((265:0.01557,200:0.00517)0.674:0.00014,((204:0.00015,(339:0.01613,(322:0.01633,268:0.01643)0.569:0.0107)0.885:0.00016)0.840:0.00527,((((((((280:0.02348,(395:0.00015,(48:0.03014,((30:0.02665,316:0.01921)0.813:0.01152,215:0.0242)0.850:0.01191)0.320:0.00016)0.912:0.02431)0.694:0.01482,(115:0.01526,364:0.08211)0.879:0.03637)0.677:0.03567,((((((162:0.06933,59:0.02113)0.991:0.08563,(308:0.02061,43:0.03488)0.894:0.04949)0.911:0.05006,(((344:0.00015,(146:0.00015,377:0.01634)0.924:0.0108)0.918:0.01069,((201:0.011,240:0.04792)1.000:0.00015,(61:0.00015,96:0.00523)0.781:0.00514)0.828:0.01056)0.809:0.00016,196:0.04505)0.213:0.00014)0.650:0.00529,(((161:0.01191,(390:0.04307,37:0.03893)0.933:0.03396)0.814:0.01401,68:0.04946)0.953:0.03303,((341:0.01127,393:0.02765)0.941:0.02238,(82:0.01112,(350:0.01141,(156:0.01636,356:0.00015)0.863:0.02214)0.946:0.02475)0.748:0.00565)0.761:0.00968)0.748:0.00836)0.927:0.0224,271:0.05902)0.753:0.00511,(((((((217:0.03796,379:0.00016)0.973:0.05805,(299:0.08963,(382:0.06426,((317:0.00016,((205:0.00532,264:0.03867)0.939:0.01605,(194:0.03374,(32:0.01052,(348:0.02212,157:0.02743)1.000:0.00014)0.868:0.02793)0.745:0.00531)0.336:0.01061)0.789:0.00604,334:0.02104)0.598:0.01527)0.687:0.00354)0.836:0.01564)0.811:0.01617,(292:0.06237,84:0.02159)0.934:0.04776)0.864:0.02103,301:0.06716)0.698:0.0046,272:0.00539)0.809:0.0115,88:0.05965)0.860:0.01208,(276:0.01065,279:0.03443)0.891:0.01124)0.090:0.00014)0.924:0.03938)0.953:0.05227,281:0.02828)0.691:0.00622,25:0.01213)0.727:0.00397,((261:0.01613,((147:0.01555,20:0.00016)0.967:0.02125,(107:0.01089,349:0.03426)0.757:0.00478)0.750:0.00518)0.799:0.0052,(259:0.01616,63:0.01053)0.764:0.00523)0.792:0.00511)1.000:0.00016,(72:0.05949,(1:0.01425,67:0.0377)0.751:0.00762)0.867:0.01609)0.807:0.00507,((49:0.01645,116:0.01633)0.736:0.00514,(398:0.00515,(((180:0.04458,99:0.0328)0.913:0.02521,(((410:0.05589,(((150:0.04425,(170:0.03163,((250:0.00693,331:0.00435)1.000:0.10845,357:0.01319)0.850:0.0225)0.879:0.02887)0.749:0.00795,(((((23:0.00919,248:0.08024)0.405:0.03691,(358:0.05635,369:0.07223)0.978:0.09469)0.888:0.05975,(234:0.07249,8:0.00016)0.712:0.01829)0.976:0.07916,(((275:0.094,(((114:0.0269,302:0.02202)0.985:0.06964,(213:0.06889,42:0.03436)0.415:0.01928)0.795:0.02064,((110:0.05188,342:0.01457)0.967:0.08524,((123:0.02756,343:0.0481)0.800:0.01738,((298:0.03283,(124:0.02507,6:0.03351)0.781:0.01076)0.939:0.03194,309:0.04124)0.820:0.01321)0.985:0.0961)0.928:0.06559)0.902:0.03886)0.684:0.03217,373:0.06838)0.909:0.03592,((290:0.02673,380:0.00015)1.000:0.16099,(((90:0.09952,192:0.10171)0.679:0.01316,(326:0.03972,45:0.09053)0.965:0.05309)0.115:0.00014,(375:0.00015,(221:0.00071,278:0.05255)1.000:0.08313)1.000:0.10921)0.623:0.0222)0.892:0.03509)0.465:0.00015)0.980:0.05443,(((306:0.08813,385:0.14214)0.269:0.00862,((256:0.01776,(273:0.07543,69:0.01333)0.591:0.02343)0.883:0.02549,((132:0.02365,219:0.01597)0.897:0.02388,(100:0.01243,50:0.0237)0.226:0.01766)0.961:0.04348)0.848:0.01577)0.998:0.08323,(241:0.23207,(130:0.24778,(53:0.12887,(129:0.07692,318:0.01288)0.900:0.04845)0.817:0.02143)0.888:0.05464)0.657:0.01537)0.822:0.01876)0.828:0.01549)0.773:0.01019,((98:0.12681,((148:0.0294,391:0.00571)0.989:0.07803,(389:0.10107,(252:0.00014,362:0.01104)0.964:0.06682)0.834:0.03217)0.762:0.0152)0.524:0.0181,(0:0.0483,(135:0.01151,(300:0.0175,(274:0.04561,((((166:0.02935,355:0.00015)0.833:0.00565,41:0.00014)0.807:0.00586,(226:0.01038,92:0.0044)0.792:0.00425)0.961:0.03236,((360:0.01752,7:0.0182)0.748:0.00495,(368:0.02316,288:0.01783)0.759:0.00622)0.707:0.00573)0.841:0.00015)0.949:0.02275)0.745:0.00559)0.855:0.02344)0.876:0.03532)0.885:0.02567)0.752:0.00645)0.782:0.00969,(((((((178:0.01576,(230:0.02704,64:0.02146)0.869:0.0108)0.809:0.01014,((122:0.00448,354:0.0225)0.855:0.01127,(333:0.01086,406:0.01648)0.748:0.00433)0.789:0.00624)0.171:0.00516,((416:0.04298,(400:0.01045,74:0.01051)0.923:0.00014)0.862:0.02166,(307:0.04097,(260:0.03574,335:0.0434)0.747:0.00875)0.916:0.02837)0.843:0.00987)0.804:0.00016,((237:0.09447,((370:0.01631,(319:0.04803,(60:0.01986,405:0.01742)0.560:0.01574)0.898:0.01971)0.918:0.01584,(384:0.02116,(245:0.01047,(177:0.0051,(183:0.03226,413:0.00014)0.826:0.00518)0.777:0.00501)0.923:0.0158)0.622:0.00016)0.685:0.00099)0.224:0.02406,((22:0.03142,5:0.06696)0.870:0.03448,47:0.0347)0.763:0.01052)0.847:0.01209)0.743:0.00534,((((62:0.00137,(121:0.00016,78:0.04376)1.000:0.10609)0.942:0.0378,(311:0.05626,407:0.06902)0.944:0.04614)0.703:0.00608,(((188:0.01993,202:0.02611)0.914:0.02118,(328:0.0273,337:0.00015)0.815:0.01019)0.852:0.01169,(330:0.03441,((386:0.13035,(392:0.00544,(321:0.02191,4:0.01061)0.763:0.0052)0.932:0.00014)0.671:0.01096,145:0.01556)0.829:0.01073)0.735:0.00529)0.840:0.01052)0.849:0.01531,(262:0.0683,((310:0.05551,((83:0.01296,(127:0.01909,212:0.01393)0.090:0.00499)0.876:0.01352,(104:0.00014,171:0.01061)0.895:0.01717)0.877:0.02683)0.940:0.03929,(119:0.0152,179:0.00197)0.889:0.02843)0.066:0.01551)0.839:0.01374)0.820:0.01069)0.869:0.01061,(((293:0.01741,168:0.04514)0.046:0.01491,345:0.03334)0.248:0.01629,(31:0.04727,97:0.04999)0.915:0.03556)0.811:0.01631)0.010:0.00016,(((94:0.0671,(108:0.00014,229:0.06991)0.630:0.01827)0.982:0.06031,(143:0.02201,((((((198:0.02745,(140:0.14724,75:0.02831)0.817:0.0209)0.851:0.01902,(((282:0.06783,54:0.00015)0.952:0.03641,((313:0.03746,80:0.00524)0.872:0.0215,2:0.07468)0.862:0.02589)0.916:0.03914,((367:0.0099,(((128:0.0425,((111:0.06727,11:0.00495)0.974:0.02953,283:0.02606)0.504:0.00357)0.862:0.02044,(289:0.04546,(399:0.00319,((((152:0.00014,19:0.06307)0.992:0.03752,154:0.00016)0.786:0.00014,134:0.06945)0.997:0.06109,51:0.00014)0.994:0.04556)0.353:0.00583)0.482:0.00828)0.933:0.03536,112:0.07957)0.734:0.00733)0.962:0.08492,403:0.10375)0.869:0.0525)0.894:0.03949)0.645:0.00925,((((287:0.00534,15:0.05518)0.920:0.03189,(((304:0.00508,409:0.00015)0.991:0.00014,(120:0.00015,(57:0.04309,56:0.0156)0.759:0.00015)0.902:0.01019)0.339:0.01644,173:0.094)0.787:0.01131)1.000:0.07731,(236:0.00625,((26:0.04569,(((351:0.005,(27:0.03624,(137:0.01569,(314:0.00015,408:0.03277)0.991:0.03257)0.806:0.00498)0.851:0.00588)0.928:0.01791,((133:0.04374,(227:0.00527,(412:0.00014,(175:0.00507,((95:0.01566,210:0.00014)0.438:0.01045,191:0.00016)0.781:0.00518)0.815:0.00508)0.859:0.01021)0.745:0.00667)0.735:0.01956,((((12:0.01588,415:0.01701)0.121:0.03139,(73:0.04886,(17:0.00016,(46:0.02083,378:0.01021)0.886:0.01027)0.785:0.019)0.719:0.02118)0.774:0.01959,329:0.01522)0.777:0.01121,(((286:0.00722,(394:0.01596,(372:0.00015,225:0.0446)0.884:0.0109)0.929:0.02558)0.584:0.00985,218:0.02283)0.888:0.01478,159:0.02121)0.739:0.00866)0.851:0.01129)0.728:0.00602)0.866:0.01998,93:0.04869)0.604:0.00297)0.648:0.01633,199:0.06704)0.788:0.01956)0.371:0.01052)0.827:0.01491,((244:0.0262,(126:0.00015,163:0.03192)0.984:0.04377)0.817:0.01306,((216:0.00014,(86:0.02257,(21:0.01127,34:0.01066)0.859:0.01088)0.622:0.017)0.998:0.19434,(233:0.00244,(182:0.01898,(239:0.02877,267:0.00015)0.839:0.01438)0.999:0.09419)0.975:0.15234)0.877:0.07457)0.893:0.0244)0.821:0.02013)0.998:0.10422,(195:0.10508,((249:0.0368,(336:0.04596,((263:0.02407,(277:0.01295,190:0.03788)0.823:0.01671)0.698:0.0068,197:0.01756)0.309:0.01631)0.860:0.01866)0.926:0.02656,(303:0.04293,(113:0.04423,347:0.04295)0.930:0.03972)0.885:0.02484)0.701:0.00015)0.902:0.03629)0.841:0.02905,(246:0.00014,(125:0.03009,184:0.0229)0.998:0.07478)0.999:0.10301)0.936:0.04978,((247:0.04204,((((((238:0.01393,(109:0.01081,39:0.02762)0.769:0.00519)0.758:0.00702,(257:0.01539,85:0.07408)0.746:0.00558)0.755:0.01039,(363:0.04294,155:0.00015)0.943:0.02426)0.894:0.01745,266:0.00586)0.948:0.03346,55:0.02705)0.739:0.00453,203:0.00015)0.855:0.02077)0.995:0.07638,327:0.00745)0.921:0.03692)0.553:0.01549)0.970:0.05544)0.858:0.02855,338:0.08163)0.892:0.03304)0.759:0.00673)0.945:0.02495,((((((((102:0.04317,36:0.02415)0.964:0.03758,65:0.00505)0.822:0.01078,366:0.00016)0.811:0.01537,(315:0.01071,((151:0.0,160:0.0):0.00016,340:0.00014)0.842:0.01037)0.951:0.02712)0.724:0.00057,(185:0.04527,(207:0.01304,76:0.00341)0.949:0.03474)0.845:0.0196)0.871:0.0106,(371:0.02805,(164:0.0104,242:0.02179)0.758:0.0052)0.771:0.00538)0.841:0.01097,174:0.13953)0.831:0.01033,(144:0.01866,(3:0.01578,312:0.00015)0.785:0.00532)0.780:0.00615)0.752:0.00572)0.667:0.00244)0.268:0.00339,((101:0.04199,77:0.00334)0.965:0.0345,((14:0.01106,294:0.00502)0.891:0.01811,(285:0.01062,397:0.01076)0.758:0.00896)0.163:0.01034)0.850:0.01331)0.563:0.00537)0.800:0.00519)0.930:0.00016)0.759:0.01023)1.000:0.00014)0.850:0.00015,(243:0.03373,220:0.01032)0.888:0.011)0.540:0.00014,(189:0.02629,(((139:0.0155,186:0.01757)0.420:0.01444,(((((((165:0.0059,58:0.03297)0.779:0.02132,((222:0.01678,(323:0.02243,44:0.04081)0.819:0.01102)0.063:0.00015,(106:0.03989,149:0.02047)0.775:0.01298)0.706:0.0074)0.957:0.03281,((((258:0.04247,87:0.0123)0.500:0.01067,235:0.00735)0.645:0.00296,208:0.00505)1.000:0.00015,((18:0.00454,(((10:0.04233,(414:0.00016,(142:0.01127,66:0.03479)0.756:0.00498)0.726:0.00685)0.486:0.01639,181:0.00014)0.784:0.00501,(167:0.01463,(320:0.00885,402:0.00881)0.791:0.00014)0.839:0.01499)0.773:0.00524)0.893:0.01079,(169:0.00517,(295:0.01586,297:0.03792)0.262:0.00016)0.778:0.00521)0.818:0.00528)0.764:0.01062)0.767:0.00486,70:0.00512)0.766:0.00495,(((332:0.00016,((325:0.01591,(383:0.00014,(361:0.01642,(138:0.04133,(158:0.0036,224:0.00657)0.840:0.01972)0.769:0.00881)0.777:0.00496)0.882:0.01036)0.752:0.00492,(24:0.03974,((((254:0.00541,(251:0.00015,(324:0.02187,((117:0.0052,(374:0.03165,270:0.02362)0.731:0.00708)0.791:0.00525,13:0.01621)0.757:0.00511)0.607:0.01283)0.889:0.0192)0.852:0.01583,305:0.01647)0.948:0.00015,211:0.00015)0.419:0.00016,(103:0.01686,209:0.05269)0.861:0.01595)0.937:0.01635)0.756:0.00523)0.878:0.01048)0.776:0.00238,(365:0.03251,((38:0.04434,79:0.00014)0.758:0.00016,(296:0.043,9:0.00518)0.693:0.0162)0.508:0.00805)0.766:0.00767)0.764:0.00313,(((359:0.02181,(16:0.04469,(232:0.01621,(118:0.03421,(29:0.01612,353:0.01494)0.293:0.01034)0.864:0.01326)0.747:0.01394)0.724:0.0072)0.911:0.01681,387:0.02755)0.761:0.00523,(346:0.01957,(376:0.04072,71:0.0547)0.829:0.0181)0.750:0.00673)0.823:0.01037)0.774:0.0054)0.789:0.005,(((228:0.00529,((401:0.02214,((187:0.00532,411:0.00526)0.801:0.00583,((89:0.027,193:0.00014)0.787:0.00524,91:0.01618)0.743:0.0045)0.548:0.00548)0.825:0.016,40:0.02807)0.778:0.00992)0.824:0.01011,255:0.05012)0.966:0.00014,(352:0.01585,396:0.00014)0.784:0.02134)0.880:0.0107)0.901:0.0194,(35:0.0209,(206:0.00836,291:0.06414)0.439:0.00793)0.753:0.00846)0.763:0.00968)0.942:0.02851,28:0.0208)0.742:0.01057)0.781:0.00811)0.802:0.02029)0.750:0.01578);"""

tutorial_otu_table = """#Full OTU Counts
#OTU ID	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636	Consensus Lineage
0	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
1	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
2	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
3	2	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
4	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
5	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
6	0	0	0	0	0	0	0	1	0	Root;Bacteria;Actinobacteria;Actinobacteria
7	0	0	2	0	0	0	0	0	2	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
8	1	1	0	2	4	0	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Lactobacillaceae;Lactobacillus
9	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
10	0	1	0	0	0	0	0	0	0	Root;Bacteria
11	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
12	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes
13	1	0	0	1	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
14	0	0	1	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
15	0	0	0	0	1	0	0	0	0	Root;Bacteria
16	1	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
17	0	0	0	1	0	0	4	10	37	Root;Bacteria;Bacteroidetes
18	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
19	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
20	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
21	0	0	0	0	0	0	2	3	2	Root;Bacteria;Bacteroidetes
22	0	0	0	0	2	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
23	14	1	14	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Lactobacillaceae;Lactobacillus
24	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
25	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
26	0	0	0	0	0	0	0	1	1	Root;Bacteria;Bacteroidetes
27	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
28	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
29	6	0	4	0	2	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
30	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
31	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
32	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
33	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
34	0	0	0	0	0	0	8	10	2	Root;Bacteria
35	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
36	1	0	1	0	0	0	0	1	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
37	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
38	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
39	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
40	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
41	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
42	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
43	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
44	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
45	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Erysipelotrichi";"Erysipelotrichales";Erysipelotrichaceae;Coprobacillus
46	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
47	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
48	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
49	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
50	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
51	0	1	0	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
52	0	2	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
53	0	0	0	0	0	0	2	0	1	Root;Bacteria;Proteobacteria;Deltaproteobacteria
54	0	0	0	0	0	0	5	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
55	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
56	0	0	0	0	0	1	0	0	0	Root;Bacteria;Bacteroidetes
57	0	0	0	0	0	0	0	1	0	Root;Bacteria
58	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
59	0	0	0	0	0	0	0	0	1	Root;Bacteria;Deferribacteres;Deferribacteres;Deferribacterales;Deferribacteraceae;Mucispirillum
60	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
61	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
62	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
63	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
64	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
65	0	0	0	6	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
66	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
67	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
68	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
69	0	0	1	0	0	0	0	0	0	Root;Bacteria
70	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
71	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
72	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
73	0	0	0	0	0	5	0	0	0	Root;Bacteria;Bacteroidetes
74	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
75	1	0	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
76	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
77	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
78	1	0	1	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
79	2	3	8	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
80	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
81	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
82	0	0	0	0	0	2	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
83	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
84	1	0	0	0	0	0	0	2	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae";Ruminococcus
85	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
86	0	0	0	0	0	0	0	1	0	Root;Bacteria
87	0	0	1	0	0	2	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
88	0	0	0	0	0	0	0	1	0	Root;Bacteria
89	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
90	0	0	0	9	0	0	3	0	0	Root;Bacteria;Firmicutes;"Erysipelotrichi";"Erysipelotrichales";Erysipelotrichaceae;Turicibacter
91	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";Butyrivibrio
92	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
93	0	0	0	0	0	0	2	1	0	Root;Bacteria;Bacteroidetes
94	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
95	0	0	0	2	0	0	0	0	0	Root;Bacteria;Bacteroidetes
96	0	0	0	1	0	1	0	1	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
97	0	0	0	0	0	1	0	0	0	Root;Bacteria
98	0	0	0	0	0	0	0	1	0	Root;Bacteria
99	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
100	0	0	0	1	0	0	0	0	0	Root;Bacteria
101	0	0	0	3	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
102	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
103	0	1	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
104	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
105	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
106	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
107	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
108	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;Incertae Sedis XIII;Anaerovorax
109	0	0	0	1	0	0	1	5	2	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
110	0	0	0	0	0	2	0	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae;Olsenella
111	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
112	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
113	0	0	0	0	0	1	0	0	0	Root;Bacteria
114	0	0	0	0	0	1	0	0	0	Root;Bacteria
115	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
116	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
117	1	0	2	0	0	6	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
118	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
119	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
120	1	3	1	2	1	9	2	4	5	Root;Bacteria;Bacteroidetes
121	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
122	0	0	0	1	0	2	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
123	0	0	0	0	0	0	1	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
124	0	0	0	0	0	0	1	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
125	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes
126	0	0	2	0	0	0	0	1	0	Root;Bacteria
127	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
128	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
129	0	0	0	1	0	0	0	0	0	Root;Bacteria
130	0	0	0	0	5	2	0	0	0	Root;Bacteria;Proteobacteria;Epsilonproteobacteria;Campylobacterales;Helicobacteraceae;Helicobacter
131	0	0	1	3	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
132	0	0	0	0	1	0	0	0	0	Root;Bacteria
133	0	0	1	0	0	0	0	0	0	Root;Bacteria
134	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
135	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
136	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
137	0	0	0	0	0	0	0	1	0	Root;Bacteria
138	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
139	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
140	0	0	0	0	0	0	1	3	0	Root;Bacteria
141	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
142	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
143	0	0	1	0	0	0	0	0	0	Root;Bacteria
144	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
145	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
146	1	0	0	0	2	0	2	0	3	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
147	0	1	0	1	1	0	0	0	3	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
148	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
149	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
150	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
151	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
152	0	0	0	1	0	0	1	2	19	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
153	0	2	1	2	0	0	1	1	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
154	2	18	0	1	0	0	21	4	4	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
155	0	0	0	0	0	5	9	5	3	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
156	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
157	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
158	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
159	0	0	0	0	0	0	0	1	1	Root;Bacteria;Bacteroidetes
160	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
161	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
162	0	0	0	0	0	3	5	2	6	Root;Bacteria;Deferribacteres;Deferribacteres;Deferribacterales;Deferribacteraceae;Mucispirillum
163	0	0	0	0	0	0	0	0	1	Root;Bacteria
164	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
165	2	1	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
166	0	0	0	0	0	0	0	1	0	Root;Bacteria
167	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
168	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
169	0	2	0	7	0	0	0	2	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
170	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
171	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
172	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
173	0	0	0	0	0	1	0	0	0	Root;Bacteria
174	1	0	0	0	10	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Peptostreptococcaceae";"Peptostreptococcaceae Incertae Sedis"
175	0	0	0	0	1	0	0	0	0	Root;Bacteria;Bacteroidetes
176	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
177	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia"
178	0	0	0	2	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
179	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
180	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
181	1	4	2	6	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
182	0	0	0	0	0	1	0	0	0	Root;Bacteria
183	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia"
184	0	0	0	1	0	0	3	1	0	Root;Bacteria;Bacteroidetes
185	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
186	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
187	0	1	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
188	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
189	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
190	0	0	0	0	0	0	0	1	0	Root;Bacteria
191	2	1	10	2	24	0	0	1	1	Root;Bacteria
192	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Streptococcaceae;Streptococcus
193	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";Butyrivibrio
194	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae";Acetanaerobacterium
195	0	0	0	0	0	1	0	0	0	Root;Bacteria
196	0	0	0	0	0	1	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
197	0	1	0	0	0	0	0	0	0	Root;Bacteria
198	0	2	0	0	0	1	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales
199	0	0	0	0	0	1	1	0	0	Root;Bacteria
200	0	0	0	2	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
201	0	0	0	1	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
202	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
203	0	2	2	4	0	5	1	5	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
204	1	4	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
205	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
206	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
207	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
208	0	2	0	2	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
209	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
210	0	0	0	0	0	0	0	0	1	Root;Bacteria
211	1	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
212	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
213	0	0	0	0	0	0	0	2	0	Root;Bacteria;Firmicutes
214	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
215	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
216	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes
217	0	0	0	0	0	2	0	1	0	Root;Bacteria
218	0	0	0	0	9	1	0	0	0	Root;Bacteria;Bacteroidetes
219	0	0	0	0	1	0	0	0	0	Root;Bacteria
220	1	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
221	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes
222	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
223	0	0	0	0	0	0	0	2	2	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
224	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
225	0	2	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
226	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
227	0	1	2	0	9	1	1	1	3	Root;Bacteria;Bacteroidetes
228	16	0	0	0	12	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
229	0	0	0	0	0	1	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;Incertae Sedis XIII
230	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
231	0	19	2	0	2	0	3	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
232	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
233	0	0	0	0	1	0	0	0	0	Root;Bacteria;Bacteroidetes
234	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Lactobacillaceae;Lactobacillus
235	0	1	1	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
236	0	0	0	0	0	2	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales
237	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
238	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
239	0	0	0	0	0	1	0	0	0	Root;Bacteria
240	0	0	0	0	0	1	0	0	0	Root;Bacteria
241	0	0	0	0	0	0	2	0	0	Root;Bacteria;TM7;TM7_genera_incertae_sedis
242	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
243	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
244	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
245	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
246	0	0	0	0	0	0	0	1	0	Root;Bacteria
247	0	0	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
248	1	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Lactobacillaceae;Lactobacillus
249	1	0	0	0	0	0	0	0	0	Root;Bacteria
250	1	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
251	0	0	0	1	4	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
252	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
253	0	0	0	0	2	0	0	5	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
254	11	13	6	13	2	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
255	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
256	0	0	0	0	0	0	1	0	0	Root;Bacteria
257	0	0	0	0	0	0	5	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
258	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
259	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
260	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
261	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
262	0	1	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";Bryantella
263	0	0	0	0	1	0	0	0	0	Root;Bacteria
264	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
265	0	0	0	0	0	2	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
266	0	0	0	2	0	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
267	1	0	0	5	17	20	0	0	0	Root;Bacteria
268	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
269	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
270	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
271	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
272	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
273	0	0	0	0	0	0	1	0	0	Root;Bacteria
274	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
275	0	0	0	0	0	0	1	0	0	Root;Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Verrucomicrobiaceae;Akkermansia
276	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
277	1	0	0	0	0	0	0	0	0	Root;Bacteria
278	0	0	0	0	0	1	0	0	0	Root;Bacteria
279	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
280	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
281	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
282	0	0	0	0	0	0	2	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
283	0	0	0	0	0	0	2	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
284	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
285	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
286	0	2	3	1	4	0	5	0	4	Root;Bacteria;Bacteroidetes
287	0	0	0	0	0	0	1	1	1	Root;Bacteria;Bacteroidetes
288	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
289	0	0	0	0	3	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
290	0	0	0	0	0	0	0	0	2	Root;Bacteria;Firmicutes;"Bacilli";Bacillales;"Staphylococcaceae";Staphylococcus
291	0	0	0	0	1	0	0	0	0	Root;Bacteria
292	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
293	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
294	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
295	29	1	10	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
296	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
297	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
298	0	0	0	0	0	0	1	0	0	Root;Bacteria;Actinobacteria;Actinobacteria
299	0	0	0	0	0	0	1	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
300	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia"
301	0	0	0	0	0	0	2	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
302	0	0	0	0	0	1	0	0	0	Root;Bacteria
303	0	0	0	0	0	0	0	0	1	Root;Bacteria
304	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes
305	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
306	0	0	0	0	0	0	0	0	1	Root;Bacteria
307	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
308	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae";"Ruminococcaceae Incertae Sedis"
309	0	0	0	1	0	0	0	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae;Denitrobacterium
310	0	0	1	0	0	0	0	0	0	Root;Bacteria
311	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
312	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
313	0	1	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
314	0	0	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
315	1	3	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
316	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
317	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
318	0	0	0	0	0	1	0	0	0	Root;Bacteria;Proteobacteria
319	0	2	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
320	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
321	0	0	0	0	0	0	0	0	1	Root;Bacteria
322	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
323	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
324	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
325	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
326	0	0	0	0	4	0	0	0	2	Root;Bacteria;Firmicutes;"Erysipelotrichi";"Erysipelotrichales";Erysipelotrichaceae;Erysipelotrichaceae Incertae Sedis
327	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes
328	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
329	2	2	0	1	0	0	0	0	0	Root;Bacteria;Bacteroidetes
330	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes
331	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes
332	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
333	0	0	0	0	0	6	0	3	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
334	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
335	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
336	0	0	1	0	0	0	0	0	0	Root;Bacteria
337	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
338	0	0	0	0	0	0	0	1	0	Root;Bacteria
339	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
340	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
341	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
342	0	0	0	0	0	1	0	0	0	Root;Bacteria
343	0	0	0	0	0	0	0	0	1	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
344	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
345	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
346	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
347	0	0	0	1	0	0	0	0	0	Root;Bacteria
348	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
349	0	0	0	0	0	0	1	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
350	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
351	0	0	0	0	2	2	1	4	1	Root;Bacteria;Bacteroidetes
352	3	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
353	0	4	4	0	1	2	0	2	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
354	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
355	0	0	0	0	0	0	0	1	0	Root;Bacteria
356	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
357	0	0	0	4	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
358	0	0	1	0	0	0	0	0	0	Root;Bacteria
359	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
360	0	0	1	0	0	0	0	1	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
361	2	0	2	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
362	1	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
363	0	0	0	0	0	1	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae
364	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
365	0	0	0	0	0	2	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
366	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";Roseburia
367	0	0	0	0	1	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
368	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
369	0	0	0	0	0	1	0	0	0	Root;Bacteria
370	2	1	0	5	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
371	1	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
372	0	1	0	0	0	0	0	0	0	Root;Bacteria
373	0	1	0	0	0	0	3	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;Clostridiaceae;"Clostridiaceae 1";Clostridium
374	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
375	0	0	0	0	0	0	4	0	0	Root;Bacteria;Firmicutes;"Erysipelotrichi";"Erysipelotrichales";Erysipelotrichaceae;Erysipelotrichaceae Incertae Sedis
376	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
377	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
378	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
379	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
380	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Bacilli";Bacillales;"Staphylococcaceae";Staphylococcus
381	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
382	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
383	4	9	0	2	0	0	0	2	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
384	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
385	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";"Carnobacteriaceae";"Carnobacteriaceae 1"
386	0	0	1	0	0	0	0	0	0	Root;Bacteria
387	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
388	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
389	0	1	0	0	0	0	0	0	0	Root;Bacteria
390	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
391	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes
392	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
393	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
394	0	0	1	0	0	0	0	0	0	Root;Bacteria
395	1	1	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
396	2	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
397	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
398	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
399	0	0	0	0	0	0	13	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
400	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
401	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
402	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
403	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Prevotellaceae
404	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
405	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
406	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
407	1	0	0	0	0	4	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
408	1	5	3	2	0	0	0	0	1	Root;Bacteria;Bacteroidetes
409	0	0	0	0	0	0	0	1	1	Root;Bacteria;Bacteroidetes
410	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
411	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
412	0	0	0	0	2	0	0	0	0	Root;Bacteria;Bacteroidetes
413	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
414	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
415	0	0	0	0	0	7	0	2	2	Root;Bacteria;Bacteroidetes
416	0	1	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales"""

unifrac_dmtx = """	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636
PC.354	0	0.623372333874	0.637204555763	0.596124964451	0.59620370837	0.736037238732	0.790494137353	0.70551354446	0.758236466255
PC.355	0.623372333874	0	0.613486299912	0.631198829777	0.672502144191	0.766682588739	0.732473435165	0.738434734323	0.687363838411
PC.356	0.637204555763	0.613486299912	0	0.691902878696	0.742922266508	0.740684448418	0.78590671873	0.723519724117	0.754301399992
PC.481	0.596124964451	0.631198829777	0.691902878696	0	0.668731241507	0.690956885033	0.652257003318	0.652311676944	0.661223123598
PC.593	0.59620370837	0.672502144191	0.742922266508	0.668731241507	0	0.74048244192	0.760255173866	0.761741286862	0.751777287666
PC.607	0.736037238732	0.766682588739	0.740684448418	0.690956885033	0.74048244192	0	0.73093223871	0.688768250911	0.729047763792
PC.634	0.790494137353	0.732473435165	0.78590671873	0.652257003318	0.760255173866	0.73093223871	0	0.583577161959	0.596033900207
PC.635	0.70551354446	0.738434734323	0.723519724117	0.652311676944	0.761741286862	0.688768250911	0.583577161959	0	0.6271222318
PC.636	0.758236466255	0.687363838411	0.754301399992	0.661223123598	0.751777287666	0.729047763792	0.596033900207	0.6271222318	0"""


big_tree_trimmed = """((((((((170836:0.10013,170502:0.07273)0.925.12:0.13018,96944:0.1819)0.937.11:0.01898,(((413456:0.04077,550248:0.04099)0.265.2:0.08525,(((115098:0.01154,274084:0.07672)1.000.397:0.11673,178347:0.08887)0.993.26:0.01497,(198054:0.06647,262375:0.06851)0.957.28:0.00965)0.607.2:0.10815)0.685.2:0.1377,470494:0.30841)0.834.27:0.01007)0.811.20:0.02003,((170335:0.07287,366044:0.11078)0.824.14:0.11839,((((((((234036:0.04224,205981:0.02934)0.963.73:0.00778,(176600:0.03643,277259:0.01351)0.934.49:0.0347)0.159.4:0.00252,(261177:0.03638,213945:0.025)0.970.69:0.0151)0.738.20:0.00347,((((264496:0.04143,430219:0.05082)0.437.4:0.00014,((322991:0.07643,177453:0.0354)0.754.18:0.00015,((((135567:0.02031,189531:0.04196)0.952.65:0.0074,185326:0.01891)0.382.6:0.00195,178596:0.06164)0.921.52:0.01612,215756:0.03762)0.715.11:0.01947)0.785.20:0.01174)0.930.48:0.00543,204144:0.06385)0.902.32:0.00667,447694:0.04416)0.367.4:0.00225)0.791.24:0.0032,(303491:0.01932,188969:0.0365)0.985.87:0.02311)0.827.18:0.00914,(169379:0.01501,(175413:0.0037,199177:0.02669)0.813.19:0.01178)0.639.8:0.13019)0.242.6:0.1051,((((((((577170:0.03957,271150:0.04157)0.383.2:0.1141,270523:0.23989)0.791.25:0.00353,555222:0.01108)0.884.46:0.00497,469842:0.02238)0.721.9:0.00801,469834:0.02418)0.957.91:0.01978,469746:0.03884)0.913.54:0.00377,129213:0.02293)1.000.1329:0.10549,((291090:0.06289,470806:0.00976)0.934.51:0.02333,249661:0.02579)0.979.100:0.05448)0.990.115:0.02255)0.957.93:0.08197,((((104780:0.01315,335530:0.0201)0.876.42:0.00362,(230364:0.03752,232099:0.01534)0.353.8:0.00553)0.073.7:0.00184,233478:0.02266)0.980.87:0.03618,569158:0.04516)0.860.40:0.10808)0.106.4:0.01914)0.947.66:0.34886)0.936.62:0.0251,((175843:0.03873,(267123:0.0162,332547:0.03943)0.699.8:0.03441)0.980.120:0.03055,461524:0.11587)0.903.64:0.2595)0.966.104:0.02159,423647:0.47617)0.940.101:0.01249,(((261590:0.1639,456006:0.12923)0.975.129:0.14243,((212758:0.1231,(227950:0.0258,215086:0.06154)0.759.21:0.0368)0.958.114:0.10557,((130890:0.02834,568692:0.0779)0.770.27:0.12789,((((170555:0.01177,(235071:0.00111,266512:0.0044)0.995.288:0.00702)0.991.207:0.0779,((((((258191:0.01004,(167034:0.0255,232326:0.00732)0.764.28:0.00156)0.882.61:0.00224,452682:0.01087)0.950.114:0.01244,(169479:0.01014,((276404:0.00938,192930:0.02166)0.923.79:0.0092,273545:0.01716)0.919.82:0.03896)0.612.12:0.02126)0.869.65:0.00458,270519:0.07696)0.981.140:0.02557,(((172136:0.01725,((355178:0.01902,(265813:0.02289,344329:0.00974)0.927.88:0.00717)0.940.106:0.00315,258899:0.02779)0.871.68:0.00259)0.899.69:0.00556,318370:0.02277)0.894.65:0.0093,232406:0.0129)0.998.419:0.01212)0.912.64:0.01427,(((261802:0.01692,134096:0.01464)0.927.90:0.01386,230534:0.00771)0.862.65:0.03192,(263106:0.01257,304013:0.04149)0.754.30:0.00517)0.672.7:0.03113)0.988.199:0.0621)0.673.11:0.01368,((322112:0.10187,((((((((((((195470:0.01109,((((131793:0.00908,(328536:0.02707,276663:0.01937)0.993.226:0.02241)0.960.129:0.00438,231146:0.01846)0.970.126:0.01052,(((229073:0.00902,(232638:0.01374,(132092:0.02011,279026:0.00734)0.861.71:0.00237)0.322.6:0.00316)0.774.26:0.00082,129526:0.02234)0.847.44:0.00813,265503:0.0377)0.787.31:0.01026)0.309.10:0.0085,189272:0.06741)0.893.70:0.00287)0.973.134:0.00652,((191772:0.01261,183362:0.03502)0.932.83:0.00869,((183428:0.01967,(14035:0.00799,195356:0.00565)0.869.69:0.01143)0.917.87:0.00576,181683:0.03326)0.831.28:0.01182)0.970.125:0.02448)0.979.152:0.00755,194780:0.08102)0.748.29:0.01744,258202:0.08075)0.876.67:0.00593,((((178173:0.02,263339:0.01728)0.896.63:0.01278,((548361:0.04709,(((339718:0.0148,274387:0.03619)0.555.11:0.01119,((319134:0.02152,(290089:0.00882,231529:0.02009)0.769.33:0.00311)0.963.109:0.00666,((232616:0.00993,(180097:0.03497,259910:0.01704)0.902.59:0.01415)0.842.55:0.00336,260756:0.02037)0.658.10:0.00521)0.964.136:0.01115)0.976.127:0.00849,((((190047:0.03774,(269815:0.02463,174366:0.02081)0.953.105:0.00891)0.963.110:0.00854,264889:0.06178)0.872.78:0.00641,(209907:0.03452,301012:0.02669)0.758.37:0.0064)0.999.657:0.01357,(((((450047:0.01828,255362:0.00816)0.286.10:0.00443,270244:0.01569)0.998.429:0.02126,135814:0.03029)0.747.36:0.00698,((174272:0.08908,(((437151:0.00188,266778:0.02055)0.997.329:0.01363,(185989:0.03725,135493:0.03486)0.914.73:0.00542)0.965.126:0.00574,312226:0.0519)0.627.10:0.00497)0.901.69:0.0039,(342948:0.03853,(338813:0.01625,180288:0.0422)0.419.8:0.02626)0.986.172:0.02407)1.000.2503:0.01858)0.207.7:0.00271,((319219:0.02396,((233411:0.00737,285281:0.02701)0.779.21:0.00893,186526:0.01256)0.977.149:0.03485)0.988.207:0.01568,173417:0.02548)0.803.31:0.01311)0.817.31:0.00751)0.902.60:0.00447)0.732.27:0.00254)0.790.23:0.00229,((302407:0.01902,553492:0.01131)0.982.181:0.0133,288931:0.0461)0.769.34:0.01051)0.747.37:0.01002)0.466.14:0.00341,(135341:0.08256,194202:0.03295)0.954.107:0.00659)0.960.127:0.008,(((135428:0.0568,(273383:0.02921,338141:0.0319)0.769.35:0.0431)0.881.70:0.00481,166592:0.12932)0.771.38:0.03068,((193626:0.04736,(346275:0.01287,268581:0.04769)0.982.186:0.0339)0.768.43:0.00183,((((274578:0.06517,403497:0.03884)0.967.131:0.00918,(164638:0.01792,308832:0.01278)0.963.113:0.00667)0.873.70:0.00492,231169:0.02881)0.969.130:0.00863,(266058:0.06419,(((412586:0.01293,271856:0.03595)0.757.47:0.00439,322062:0.03984)0.899.76:0.013,534050:0.04296)0.723.24:0.0053)0.911.82:0.00732)0.914.77:0.00504)0.985.158:0.011)0.831.27:0.00677)0.914.78:0.01057)0.892.58:0.00552,(((320141:0.06784,(194822:0.01044,326759:0.02963)0.831.29:0.0074)0.981.151:0.00871,(((132165:0.02094,(343103:0.01541,230541:0.00672)0.965.135:0.00689)0.974.139:0.00928,187274:0.01062)0.965.136:0.00861,325636:0.03123)0.203.4:0.01001)0.908.89:0.01865,343906:0.03179)0.921.87:0.00515)0.949.114:0.00497,(((328453:0.03326,(172622:0.00371,386087:0.01456)0.336.16:0.0041)0.826.26:0.00188,(309206:0.03467,271097:0.01673)1.000.2484:0.00282)0.913.84:0.01587,192659:0.03243)0.985.153:0.01882)0.990.201:0.00015,((275869:0.07902,(172705:0.01959,230573:0.01377)0.999.681:0.05317)0.252.4:0.02525,271439:0.06032)0.833.41:0.00145)0.916.82:0.00554,(423411:0.07406,(((((177802:0.01315,269872:0.01788)0.999.682:0.03959,(131203:0.02688,264373:0.01833)0.986.181:0.0214)0.835.38:0.01945,272454:0.03294)0.996.316:0.01321,232142:0.03761)0.489.4:0.00287,(((350876:0.01991,274021:0.02901)0.601.9:0.00717,315527:0.01247)0.968.129:0.01233,((262047:0.0214,130526:0.02995)0.888.62:0.00868,(332768:0.01281,296869:0.05014)0.348.6:0.02613)0.826.34:0.02606)0.971.165:0.00793)0.910.80:0.02542)0.905.83:0.02439)0.895.74:0.00369,(195445:0.04495,339886:0.0318)0.983.173:0.02211)0.387.13:0.00159,232321:0.15374)0.803.32:0.00334,(((((((322487:0.02133,316210:0.03285)0.977.168:0.02558,342570:0.06745)0.900.74:0.00592,(((261606:0.02821,337617:0.02183)0.939.112:0.00749,(260205:0.02858,231684:0.02209)0.373.6:0.00559)0.989.193:0.01987,(275847:0.03082,182470:0.02718)0.930.99:0.08341)0.990.209:0.0113)0.887.62:0.00525,((342458:0.03995,263876:0.03166)0.702.18:0.02281,411144:0.02783)0.993.237:0.01961)0.982.200:0.01175,163240:0.11626)0.974.148:0.01179,(258750:0.02047,186426:0.05109)0.998.453:0.02982)0.807.28:0.00299,241287:0.06828)0.357.12:0.0179)0.830.40:0.02751)0.978.149:0.01314,191951:0.05079)0.614.13:0.06606)0.949.121:0.01032,((((231955:0.0086,264695:0.00727)0.960.133:0.08278,((229754:0.04256,((267135:0.00339,170950:0.00918)0.989.196:0.00885,209866:0.02045)1.000.2644:0.0493)0.357.13:0.00159,(((407007:0.0453,269378:0.08698)0.994.259:0.02192,355041:0.08393)0.601.10:0.00426,(232601:0.0188,215193:0.0169)0.997.364:0.02804)0.462.9:0.00218)0.988.228:0.02646)0.949.125:0.03311,260352:0.16855)0.806.28:0.0169,259499:0.08749)0.940.123:0.02251)0.721.19:0.02446)0.848.58:0.03812)0.919.99:0.03152)0.979.158:0.05135,(((200058:0.10099,(351560:0.11914,((272882:0.11244,266771:0.10802)0.922.108:0.00452,(239571:0.02554,412034:0.00583)0.858.63:0.04591)0.779.30:0.01871)0.983.187:0.0512)0.278.8:0.00612,(150117:0.00581,164612:0.05052)0.588.9:0.0441)0.934.90:0.12745,468179:0.15926)0.476.7:0.02539)0.989.220:0.03427)0.692.9:0.24849);"""

big_dmtx = """	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636
PC.354	0	0.527047269	0.647737939	0.547217109	0.496056681	0.687188863	0.701120488	0.649585288	0.681434152
PC.355	0.527047269	0	0.541921911	0.533291513	0.579027153	0.651584427	0.649520247	0.620218939	0.587801465
PC.356	0.647737939	0.541921911	0	0.646436943	0.67217105	0.718922954	0.762200461	0.692259868	0.667052119
PC.481	0.547217109	0.533291513	0.646436943	0	0.625068445	0.631004563	0.636227528	0.64047777	0.588969244
PC.593	0.496056681	0.579027153	0.67217105	0.625068445	0	0.611612148	0.682412189	0.684436784	0.633293071
PC.607	0.687188863	0.651584427	0.718922954	0.631004563	0.611612148	0	0.632078639	0.680466604	0.620411354
PC.634	0.701120488	0.649520247	0.762200461	0.636227528	0.682412189	0.632078639	0	0.646910115	0.500737813
PC.635	0.649585288	0.620218939	0.692259868	0.64047777	0.684436784	0.680466604	0.646910115	0	0.601862605
PC.636	0.681434152	0.587801465	0.667052119	0.588969244	0.633293071	0.620411354	0.500737813	0.601862605	0"""


big_otu_table = """# QIIME v1.3.0 OTU table
#OTU ID	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636
104780	0	0	0	2	0	0	0	0	0
115098	0	0	0	5	17	20	0	0	0
129213	0	0	0	0	3	0	0	0	0
129526	0	16	1	0	2	0	3	0	0
130526	0	0	0	0	0	1	0	0	0
130890	1	1	1	0	0	0	0	0	0
131203	0	0	0	1	0	0	0	0	0
131793	0	0	0	1	0	0	0	0	0
132092	0	0	0	1	0	0	0	0	0
132165	1	0	0	0	0	0	0	0	0
134096	0	0	0	0	0	1	0	0	0
135341	0	0	0	0	0	1	0	0	0
135428	0	0	0	0	0	1	0	0	0
135493	11	12	4	13	0	0	0	0	0
135567	1	1	6	0	8	0	0	1	1
135814	0	0	2	0	0	0	0	0	0
14035	0	0	0	0	2	0	0	3	0
150117	0	0	0	0	0	0	0	0	1
163240	0	0	0	1	0	0	0	0	0
164612	0	0	0	0	0	0	0	0	1
164638	1	4	0	0	0	0	0	0	0
166592	0	0	0	0	0	0	0	1	0
167034	0	0	0	0	0	1	0	0	0
169379	1	2	0	1	1	6	1	4	4
169479	0	1	0	0	0	0	0	0	0
170335	0	0	0	0	0	0	3	1	0
170502	0	0	0	0	0	1	0	0	0
170555	0	0	0	1	0	0	0	0	0
170836	0	0	0	0	0	0	2	0	1
170950	0	0	0	0	0	1	0	0	0
172136	0	0	1	0	0	0	0	1	0
172622	0	0	0	1	0	0	0	0	0
172705	0	0	0	1	0	0	0	0	0
173417	0	0	0	1	0	0	0	0	0
174272	0	0	0	0	0	0	0	0	1
174366	1	0	0	0	0	0	0	0	0
175413	0	0	0	0	0	0	0	1	0
175843	0	0	0	0	0	0	0	0	1
176600	0	1	1	0	0	0	0	0	0
177453	0	0	0	2	0	0	0	0	0
177802	0	0	0	2	0	0	0	0	0
178173	0	0	1	0	0	0	0	0	0
178347	0	1	0	0	0	0	3	0	0
178596	1	5	3	2	0	0	0	0	1
180097	0	1	0	0	1	0	0	0	0
180288	0	0	1	0	0	0	0	0	0
181683	1	0	1	0	0	0	0	0	0
182470	0	0	0	0	0	0	0	0	1
183362	0	0	1	0	0	0	0	0	0
183428	1	0	0	0	0	0	0	0	0
185326	0	0	0	0	2	0	0	0	0
185989	1	0	1	0	2	0	0	0	0
186426	0	0	1	0	0	0	0	0	0
186526	0	1	0	0	0	0	0	0	0
187274	0	0	2	0	0	0	0	0	0
188969	0	0	0	0	0	0	1	0	0
189272	0	0	0	0	0	0	1	0	0
189531	1	0	4	2	16	0	0	0	0
190047	0	1	0	0	0	0	0	0	1
191772	0	0	0	0	0	1	0	0	0
191951	0	0	1	0	0	0	0	0	0
192659	0	0	0	1	0	0	0	0	0
192930	0	0	1	0	0	0	0	0	0
193626	0	1	0	1	1	0	0	0	3
194202	0	0	0	0	1	0	0	0	0
194780	0	1	0	0	1	0	0	0	0
194822	1	0	0	0	0	0	0	0	0
195356	0	0	0	0	0	0	0	1	0
195445	1	0	0	0	0	0	0	0	0
195470	1	0	0	0	0	0	0	0	0
198054	0	0	0	0	0	0	3	0	0
199177	0	0	0	0	0	0	0	1	1
200058	0	0	0	8	0	0	3	0	0
204144	0	0	3	0	1	0	2	0	2
205981	0	1	0	1	3	0	3	0	1
209866	0	0	0	0	0	1	0	0	0
209907	2	0	1	0	2	0	0	0	0
212758	0	0	0	0	0	1	0	0	0
213945	0	0	0	0	0	2	0	0	0
215086	0	0	0	0	0	1	0	0	0
215193	0	0	0	1	0	0	0	0	0
215756	0	0	0	0	4	1	0	0	1
227950	0	0	0	0	0	2	0	0	0
229073	0	0	1	1	0	0	0	0	0
229754	0	0	0	0	0	0	0	1	0
230364	0	0	0	0	0	5	9	5	3
230534	0	0	1	0	0	0	0	1	0
230541	1	2	1	0	0	0	0	0	0
230573	0	0	0	0	0	0	0	1	0
231146	0	0	0	1	0	0	0	0	0
231169	0	0	1	0	0	0	0	0	0
231529	29	1	10	0	0	0	0	0	0
231684	0	0	0	0	0	1	0	0	0
231955	0	0	0	0	1	0	0	0	0
232099	0	0	0	1	0	0	1	5	2
232142	0	0	0	2	0	0	0	0	0
232321	0	0	1	1	0	0	0	0	0
232326	0	0	1	0	0	0	0	0	1
232406	0	0	0	0	2	0	2	0	3
232601	0	0	0	0	0	0	0	1	0
232616	0	0	0	0	0	1	0	0	0
232638	0	0	1	0	0	0	0	0	0
233411	0	0	0	1	0	0	0	0	0
233478	0	2	2	4	0	5	1	5	0
234036	0	0	2	0	0	0	0	1	0
235071	0	0	0	1	0	0	0	0	0
239571	0	0	0	2	3	0	0	0	0
241287	0	0	0	1	0	0	0	0	0
249661	0	0	0	0	0	0	0	0	1
255362	1	4	2	5	0	0	0	0	0
258191	0	0	2	0	0	0	0	0	1
258202	0	0	1	0	0	0	0	1	1
258750	0	0	0	5	0	0	0	1	0
258899	1	0	0	0	0	0	0	0	0
259499	0	0	0	1	0	0	0	0	0
259910	0	0	1	0	0	0	0	0	0
260205	0	0	0	1	0	0	0	1	0
260352	0	0	0	0	0	0	1	0	1
260756	0	1	1	0	0	0	0	0	0
261177	0	0	0	0	0	1	1	0	0
261590	1	0	0	0	10	0	0	0	0
261606	0	0	0	0	0	0	0	1	0
261802	0	0	0	0	0	1	0	0	0
262047	0	1	0	0	0	0	0	0	0
262375	0	0	0	0	0	0	0	1	0
263106	0	0	0	0	0	1	0	0	0
263339	4	9	0	1	0	0	0	2	0
263876	0	3	3	0	1	2	0	2	1
264373	0	0	0	0	0	0	0	0	1
264496	2	2	0	1	0	0	0	0	0
264695	0	0	1	0	0	0	0	0	0
264889	0	0	1	0	0	0	0	0	0
265503	0	0	0	3	0	0	0	1	0
265813	0	0	0	1	0	1	0	0	0
266058	0	0	0	0	0	0	0	2	2
266512	1	0	0	0	0	0	0	0	0
266771	14	1	14	1	0	0	0	0	0
266778	0	0	0	0	0	0	1	0	0
267123	0	0	0	1	0	0	0	0	0
267135	0	0	0	0	0	0	0	1	0
268581	1	0	0	0	0	0	0	0	0
269378	1	0	0	0	0	0	0	2	0
269815	1	1	0	5	0	1	0	0	0
269872	0	0	0	0	1	0	0	0	0
270244	0	0	0	0	2	0	0	0	0
270519	0	0	0	0	0	0	0	1	0
270523	0	0	0	0	0	0	0	1	0
271097	0	0	0	0	1	0	0	0	0
271150	0	0	0	1	0	0	1	2	19
271439	0	0	0	0	1	0	0	0	0
271856	0	0	0	0	0	1	0	0	0
272454	0	0	0	1	0	0	0	0	0
272882	1	0	0	1	0	0	0	0	0
273383	0	0	0	1	0	2	0	0	0
273545	0	0	0	0	0	0	1	0	0
274021	0	0	1	0	0	0	0	0	0
274084	1	0	0	0	0	0	0	0	0
274387	0	0	1	0	0	0	0	0	0
274578	0	0	0	0	2	0	1	0	0
275847	0	0	1	0	0	0	0	0	0
275869	0	1	0	0	0	0	0	0	0
276404	0	0	0	0	0	0	0	1	0
276663	0	0	0	0	0	0	1	0	1
277259	0	1	0	0	0	0	0	0	0
279026	0	2	0	0	0	0	0	0	0
285281	0	0	0	0	0	1	0	0	0
288931	0	0	0	0	0	1	0	0	0
290089	0	1	0	7	0	0	0	2	0
291090	0	0	0	0	0	0	3	0	0
296869	0	0	2	0	0	0	0	0	0
301012	0	0	1	0	0	0	0	0	0
302407	0	1	0	0	0	0	0	0	0
303491	0	1	0	0	0	1	0	0	0
304013	0	0	1	0	0	0	0	0	0
308832	0	0	0	1	0	0	0	0	0
309206	0	0	0	0	0	0	0	1	0
312226	1	0	1	0	0	0	0	0	0
315527	0	0	0	0	0	6	0	3	0
316210	0	0	0	2	0	0	0	0	0
318370	0	0	0	1	0	1	0	0	0
319134	3	0	0	0	3	0	0	0	0
319219	0	0	0	1	0	0	0	0	0
320141	0	0	0	0	0	0	0	1	0
322062	0	0	1	0	0	0	0	0	0
322112	0	0	1	1	0	0	1	0	1
322487	0	0	0	1	0	0	0	0	0
322991	0	0	0	0	0	0	0	1	1
325636	0	0	0	1	0	0	0	1	0
326759	0	0	0	0	1	0	0	0	0
328453	0	0	0	0	0	0	1	0	0
328536	13	0	0	0	7	0	0	0	0
332547	0	0	0	0	0	0	1	0	0
332768	0	0	1	0	0	0	0	0	0
335530	0	0	0	0	0	0	5	0	0
337617	0	0	0	0	0	0	0	1	0
338141	1	0	0	0	0	0	0	0	0
338813	0	0	0	1	0	0	0	0	0
339718	1	0	0	2	0	0	0	0	0
339886	0	0	0	0	0	0	1	0	0
342458	0	0	0	0	1	0	0	0	0
342570	0	2	1	0	0	0	0	0	0
342948	1	0	0	0	0	1	0	0	0
343103	0	1	0	0	0	0	0	0	0
343906	0	0	0	0	0	0	1	0	0
344329	0	0	0	0	0	1	0	1	1
346275	1	1	0	0	0	0	0	0	0
350876	0	0	0	1	0	0	0	0	0
351560	0	0	0	0	0	1	0	0	0
355041	0	0	0	0	0	0	0	1	0
355178	0	0	1	0	0	0	0	0	0
366044	0	0	0	0	0	0	1	3	0
386087	1	0	0	0	0	0	0	0	0
403497	1	0	0	0	0	0	0	0	0
407007	0	0	0	0	0	0	0	1	0
411144	0	0	0	1	0	0	0	1	0
412034	1	0	0	0	1	0	0	0	0
412586	0	1	0	0	0	0	0	0	0
413456	1	0	0	0	0	0	0	0	0
423411	0	0	0	0	0	0	0	1	0
423647	0	0	0	0	0	0	2	0	0
430219	0	0	0	0	0	5	0	0	0
437151	0	0	0	0	0	6	0	0	0
447694	0	0	0	0	9	1	0	0	0
450047	2	0	0	0	0	0	0	0	0
452682	0	0	0	0	0	0	0	1	1
456006	0	0	0	0	0	1	1	0	0
461524	0	0	0	0	0	2	0	0	0
468179	1	0	0	0	1	0	0	1	0
469746	0	0	0	0	0	0	1	0	0
469834	0	0	0	0	0	0	1	0	0
469842	0	0	0	0	0	0	13	0	0
470494	0	0	0	0	5	2	0	0	0
470806	0	0	0	0	0	0	2	0	0
534050	1	0	0	0	0	0	0	0	0
548361	4	0	3	0	0	0	0	0	0
550248	0	0	0	0	4	0	0	0	2
553492	0	1	4	0	0	0	0	0	0
555222	0	0	0	0	0	0	1	1	0
568692	0	0	0	0	0	0	0	1	0
569158	0	0	0	0	0	0	0	0	1
577170	2	18	0	1	0	0	21	4	4
96944	0	0	0	0	0	3	5	2	6
None1	0	0	0	0	0	0	0	0	1
None10	0	0	0	0	0	1	0	0	0
None100	0	0	1	0	0	0	0	0	0
None101	0	0	0	0	0	1	0	0	0
None102	0	0	0	1	0	0	0	0	0
None103	0	0	0	0	0	1	0	0	0
None104	0	0	0	0	0	0	8	10	2
None105	0	0	1	0	0	0	0	0	0
None106	0	0	0	0	0	0	0	1	0
None107	0	0	1	0	0	0	0	0	0
None108	0	0	0	0	0	0	0	0	1
None109	0	0	0	0	0	0	0	0	1
None11	0	0	0	0	0	0	0	1	0
None110	0	0	0	0	1	0	0	0	0
None111	0	0	0	0	0	0	1	0	0
None112	0	0	0	0	0	1	0	0	0
None113	0	0	1	0	0	0	0	0	0
None114	0	0	0	1	4	0	0	0	0
None115	1	0	0	0	0	4	0	0	0
None116	0	1	0	0	0	0	0	0	0
None117	2	0	2	1	0	0	0	0	0
None118	0	0	0	4	0	0	0	0	0
None119	0	0	0	0	0	0	0	0	1
None12	0	0	0	0	0	0	1	0	0
None120	0	0	0	0	0	1	0	0	0
None121	0	0	0	0	0	1	0	0	0
None122	0	1	0	0	0	0	0	0	0
None123	1	0	0	0	0	0	0	0	0
None124	1	0	0	0	0	0	0	0	0
None125	0	0	0	1	0	0	0	0	0
None126	0	0	0	0	0	0	0	1	0
None127	0	0	0	0	0	0	0	0	1
None128	0	0	0	0	0	0	0	1	0
None129	0	0	0	0	0	0	0	1	0
None13	0	0	1	0	0	0	0	0	0
None130	0	0	0	0	0	0	0	1	0
None131	1	0	0	0	0	0	0	0	0
None132	0	0	1	0	0	0	0	0	0
None133	0	0	0	0	0	0	0	1	0
None134	1	0	0	0	0	0	0	0	0
None135	0	0	0	0	0	7	0	2	2
None136	0	0	0	0	0	0	0	0	1
None137	0	0	0	0	0	0	1	0	0
None138	0	0	0	0	0	1	0	1	0
None139	2	3	4	0	1	0	0	0	0
None14	0	0	0	0	0	1	0	0	0
None140	0	0	1	0	0	0	0	0	0
None141	0	0	0	1	0	0	0	0	0
None142	0	0	1	0	0	0	0	0	0
None143	0	0	0	1	0	0	0	0	0
None144	0	0	1	0	0	0	0	0	0
None145	1	0	0	0	0	0	0	0	0
None146	0	0	1	0	0	0	0	0	0
None147	0	0	0	0	0	0	0	0	1
None148	0	0	1	0	0	0	0	0	0
None149	0	0	0	0	0	0	0	0	1
None15	0	0	0	0	0	1	0	0	0
None150	0	0	0	1	0	0	0	0	0
None151	0	0	0	0	1	0	0	0	0
None152	0	0	0	1	0	0	0	0	0
None153	0	2	0	2	0	0	0	1	0
None154	0	0	0	0	1	0	0	0	0
None155	0	0	0	0	1	0	0	0	0
None156	0	0	0	0	1	0	0	0	0
None157	0	1	0	0	0	0	1	0	0
None158	0	0	0	0	0	0	1	0	0
None159	0	1	0	0	0	0	0	0	0
None16	0	0	0	0	0	1	0	0	0
None160	0	0	1	1	0	0	1	0	1
None161	0	0	0	0	0	1	0	0	0
None162	0	0	0	0	0	0	0	1	0
None163	0	0	0	0	0	0	0	1	0
None164	0	1	0	0	0	0	0	0	0
None165	0	0	0	0	0	0	0	0	1
None166	0	0	0	0	0	1	0	0	0
None167	0	0	0	0	0	1	0	0	0
None168	0	0	0	0	0	0	0	1	0
None169	0	1	0	0	0	0	0	0	0
None17	0	0	0	1	0	0	0	0	0
None170	0	1	0	0	0	0	0	0	1
None171	0	1	2	0	5	0	0	0	2
None172	0	0	0	1	0	0	4	10	37
None173	0	0	0	0	0	1	0	0	0
None174	0	1	0	0	0	3	0	0	0
None175	0	1	0	0	0	0	0	0	0
None176	0	0	0	0	0	1	0	0	0
None177	0	0	0	0	0	0	0	1	0
None178	0	0	0	0	0	1	0	0	0
None179	0	0	0	0	0	0	0	1	1
None18	0	0	0	0	1	0	0	0	0
None180	0	0	0	0	0	0	1	0	0
None181	0	0	0	0	0	0	1	0	0
None182	0	0	0	0	0	0	1	0	0
None183	0	0	1	0	0	0	0	0	0
None184	0	0	0	0	0	1	0	0	0
None185	0	0	0	0	0	1	0	0	0
None186	0	0	0	0	1	0	0	0	0
None187	0	0	0	0	0	0	0	0	1
None188	1	0	0	0	0	0	0	0	0
None189	0	0	2	0	0	0	0	0	0
None19	1	0	0	0	0	0	0	0	0
None190	0	1	0	0	0	0	0	0	0
None191	0	1	0	0	0	0	0	0	0
None192	0	0	0	0	0	0	0	0	1
None193	0	0	0	0	0	1	0	0	0
None194	0	0	1	0	0	0	0	0	0
None195	0	0	0	0	0	1	0	0	0
None196	3	0	0	0	0	0	0	0	0
None197	0	0	0	0	0	0	0	1	0
None198	0	0	0	0	0	0	0	0	1
None199	0	0	0	0	0	0	0	1	0
None2	1	0	0	0	0	0	0	0	0
None20	0	0	0	0	0	0	0	0	1
None200	0	0	0	0	0	0	0	1	0
None201	0	0	0	0	0	0	2	1	0
None202	0	0	0	1	0	0	0	0	0
None203	0	0	0	1	0	0	0	0	0
None204	0	0	0	0	0	1	0	0	0
None205	0	0	0	1	0	0	0	0	0
None206	0	0	1	0	0	0	0	0	0
None207	0	0	0	0	0	0	0	1	0
None208	0	0	0	0	0	0	1	1	0
None209	0	0	0	0	0	0	0	1	0
None21	0	0	0	0	0	1	0	0	0
None210	2	1	2	0	0	0	0	0	0
None211	1	1	0	0	0	0	0	0	0
None212	0	0	0	0	0	1	0	0	0
None213	0	0	0	0	0	0	0	1	1
None214	0	0	0	0	0	0	0	1	0
None215	0	1	0	0	0	0	0	0	0
None216	0	1	1	0	0	0	0	0	0
None217	0	0	0	0	1	0	0	0	0
None218	0	0	1	0	0	0	0	0	0
None219	1	0	1	0	0	0	0	0	0
None22	0	0	0	0	0	1	0	0	0
None220	0	0	0	0	0	0	2	0	0
None221	0	0	0	0	0	0	1	0	0
None222	0	0	1	0	0	0	0	0	0
None223	0	0	0	0	0	0	0	1	0
None224	0	0	0	0	0	0	0	1	0
None225	0	0	1	0	0	0	0	0	0
None226	0	0	1	0	0	0	0	0	0
None227	0	0	0	0	0	0	1	0	0
None228	0	0	0	0	0	0	0	0	1
None229	0	1	0	0	0	0	0	0	0
None23	0	0	0	0	0	0	0	0	1
None230	0	1	0	0	0	0	0	0	0
None231	0	0	0	0	0	0	0	0	1
None232	0	0	0	0	0	0	0	1	0
None233	0	1	0	0	0	0	0	0	0
None234	0	0	0	1	0	0	0	0	0
None235	0	0	0	1	0	0	0	0	0
None24	0	1	0	0	0	0	0	0	0
None25	0	0	0	1	0	0	0	0	0
None26	0	0	0	1	0	0	0	0	0
None27	0	1	0	0	0	0	0	0	0
None28	0	0	0	0	0	1	0	0	0
None29	0	0	0	0	2	2	1	4	1
None3	0	0	0	0	1	0	0	0	0
None30	0	0	0	0	0	0	2	3	2
None31	0	0	1	0	0	0	0	0	0
None32	0	0	0	0	0	0	1	0	0
None33	0	0	1	0	0	0	0	0	0
None34	0	0	0	0	0	0	0	1	0
None35	0	0	0	0	0	0	0	1	0
None36	0	0	0	0	0	0	0	1	0
None37	0	0	0	0	1	0	0	0	0
None38	0	0	1	0	0	0	0	0	0
None39	0	0	1	0	0	0	0	0	0
None4	0	0	0	1	0	0	0	0	0
None40	0	0	0	0	0	0	1	0	0
None41	0	0	0	0	0	1	0	0	0
None42	0	0	0	0	0	1	0	0	0
None43	0	0	0	0	0	0	1	0	0
None44	0	0	0	0	0	0	0	0	1
None45	0	1	0	0	0	0	0	0	0
None46	0	1	0	0	0	0	0	0	0
None47	0	1	0	0	0	0	0	0	0
None48	0	1	0	0	0	0	0	0	0
None49	0	0	0	0	0	0	1	0	0
None5	0	0	0	0	1	0	0	0	0
None50	0	0	0	0	0	1	0	0	0
None51	0	0	0	0	1	0	0	0	0
None52	0	0	0	0	0	0	0	0	1
None53	0	0	0	1	0	0	0	0	0
None54	0	1	0	0	0	0	0	0	0
None55	1	0	0	0	0	0	0	0	0
None56	0	1	0	0	0	0	0	0	1
None57	0	1	0	0	0	0	0	0	0
None58	1	0	0	0	0	0	0	0	0
None59	0	0	0	0	0	0	0	1	0
None6	0	1	0	0	0	0	0	0	0
None60	1	0	0	0	0	0	0	0	0
None61	0	0	0	0	0	1	0	0	0
None62	0	0	0	0	0	1	0	0	0
None63	0	0	1	0	0	0	0	0	0
None64	1	0	0	0	0	0	0	0	0
None65	0	0	0	0	0	0	1	0	0
None66	0	0	0	0	0	1	0	0	0
None67	0	0	0	0	0	1	0	0	0
None68	0	0	0	0	0	1	0	0	0
None69	0	0	0	1	0	0	0	0	0
None7	0	0	0	0	0	0	1	0	0
None70	0	0	0	0	0	0	0	2	0
None71	0	0	0	0	0	0	0	0	1
None72	0	0	0	0	0	0	0	1	0
None73	0	0	0	1	0	0	0	0	0
None74	0	1	1	0	0	0	0	0	0
None75	0	0	0	0	0	0	0	0	1
None76	0	0	0	0	0	0	1	0	0
None77	0	0	0	1	0	0	0	0	0
None78	0	2	0	0	0	0	0	0	0
None79	0	0	1	0	0	0	0	0	0
None8	0	0	0	0	0	0	1	0	0
None80	1	0	1	0	0	0	0	0	0
None81	0	0	0	0	0	0	0	1	0
None82	0	0	0	1	0	0	0	0	0
None83	1	0	2	0	0	0	0	0	0
None84	0	0	0	0	0	0	0	1	0
None85	1	0	0	0	0	0	0	0	0
None86	0	0	0	0	1	0	0	0	0
None87	0	0	0	0	1	0	0	0	0
None88	0	0	0	0	0	0	0	0	1
None89	0	0	0	0	0	0	2	0	0
None9	0	0	0	0	0	0	1	0	0
None90	0	0	0	0	0	1	0	0	0
None91	0	0	1	0	0	2	0	1	0
None92	0	0	0	1	0	0	0	0	0
None93	0	0	0	0	0	0	1	0	0
None94	0	0	0	1	0	0	0	0	0
None95	0	1	0	0	0	0	0	0	0
None96	0	0	0	0	0	0	1	0	0
None97	0	0	0	0	0	0	0	1	0
None98	0	0	1	0	0	0	0	0	0
None99	0	0	1	0	0	0	0	0	0"""



if __name__ == "__main__":
    main()