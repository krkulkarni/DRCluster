#!/bin/env python

__author__ = 'kulkarnik'
import argparse
import time, sys
import subprocess
import string
import wget, os
import collections
from helper import clean_pdb, amino_acids

def arg_parser():

    ## Parses all arguments
    ## For clean and reusuable options, store all flags in a 'setupoptions' file
    ## Example of setupoptions:
    """
    -align
    -p path/to/protein.fas
    -db path/to/pdbdatabase #FASTA file of all PDB entries
    -hmmer path/to/jackhmmer/executable
    -eval 1.0
    -iter 5

    ## Run rosetta_setup using python3 if possible.
    ## Ex: python3 rosetta_setup @setupoptions

    """
    paraParser = argparse.ArgumentParser(description='Setup Rosetta Run',
                                         fromfile_prefix_chars='@')
    def convert_arg_line_to_args(arg_line):
        for arg in arg_line.split():
            if not arg.strip():
                continue
            yield arg

    paraParser.convert_arg_line_to_args = convert_arg_line_to_args

    ## REQUIRED ARGUMENT
    paraParser.add_argument('-p', '--protein',
                            help="Protein FASTA file",
                            required=True)
    paraParser.add_argument('-dir', '--basedir',
                            help="Base directory",
                            default=os.getcwd())

    ## HMMER ARGUMENTS
    paraParser.add_argument('-align', '--align',
                            help="Use jackhmmer to find aligned sequences",
                            action="store_true")
    paraParser.add_argument('-db', '--pdbdb',
                            help="Name of PDB database (or other relevant DB)",
                            default=None)
    paraParser.add_argument('-hmmer', '--hmmerpath',
                            help="Path to hmmer executable",
                            default=None)
    paraParser.add_argument('-iter', '--iterations',
                            help="Number in iterations in HMMER search",
                            default=5)
    paraParser.add_argument('-eval', '--evalue',
                            help="E-value cutoff for HMMER",
                            default=1.0)

    ## ROSETTA ARGUMENTS
    paraParser.add_argument('-cpus', '--cpus',
                            help="Number of CPUS to use",
                            default=1)
    paraParser.add_argument('-n', '--numdecoys',
                            help="Number of decoys to use",
                            default=1)
    paraParser.add_argument('-rosdb', '--rosettadb',
                            help='Path to Rosetta database',
                            required=True)
    paraParser.add_argument('-pdbs', '--pdbtemplates',
                            help="List of PDB IDs to be used for templates",
                            default=None)
    paraParser.add_argument('-alnfile', '--alignmentfile',
                            help="List of alignments corresponding to PDBs",
                            default=None)

    args = paraParser.parse_args()

    return args






def execute_commands(command_array, wait=True):
    """ This function invokes the subprocess.Popen method to run system commands which are provided
        as a lists of lists, ie.   [ [ 'a.py',  '-f', 'file_name' ], [ 'b.py', '-j', '10' ] ]
        Commands are executed in the order that they are found in the list of lists and the function
        will wait for them to conclude and return their results as a list of lists.   The first element
        of each list will be the stdout output and the second element will be the stderr output from
        executed commands. Commands will be ran in the background if the optional function parameter 'wait'
        is set to False and the stderr and stdout results will not be returned.
    """
    command_results = []

    for command in command_array:

        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()

        command_results.append([stdout,stderr])

    return command_results


def runjackhmmer(fastapath,hmmer,pdbdb,iter,eval):

    ## This function accepts the path to the target fasta sequence and the fasta PDB database
    ## and creates a domain hits file (dom.hits),
    ## table of hits file (tbl.hits)
    ## jackhmmer alignment file (aln.out)
    ## and a hmmer log file (hmmerlog)
    ## see Documentation at http://hmmer.janelia.org/software for additional info

    fastaname = fastapath.split("/")[-1]
    protname = fastaname.split(".")[0]

    comm = [
         [ hmmer,
           '--domtblout', 'dom.hits', '--tblout', 'tbl.hits',
           '-A', 'aln.out', '-o', 'hmmerlog',
           '--domE', str(eval), '-E', str(eval), '-N', str(iter),
           fastapath, pdbdb
         ]
       ]

    print(' '.join(comm[0]))
    r = execute_commands(comm)
    print('result (list of lists) [0]->stdout [1]->stderr:')
    print(r)






class Alignment(object):

    ## This class holds alignment information
    ## about templates identified by jackhmmer search
    ## derived from 'aln.out' file

    def __init__(self,pdbid,chain,qstart,qalign,tstart,talign,score):
        self.pdbid = pdbid
        self.chain = chain
        self.qstart = qstart
        self.qalign = qalign
        self.tstart = tstart
        self.talign = talign
        self.score = score





class RunJM(object):

    def __init__(self, fasta, tbl, align):

        fastahandle = open(fasta, 'rt')
        # stores entire target fasta file
        self.fastafile = fastahandle.readlines()
        fastahandle.close()

        # identify protein name
        fastaname = fasta.split("/")[-1]
        self.protname = fastaname.split(".")[0]

        # create handles for tbl.out, aln.out
        self.tblhandle = open(tbl, 'rt')
        self.alnhandle = open(align, 'rt')

        # infodict will alignment information about each template
        # using template PDB ID as the key
        self.infodict = collections.OrderedDict()
        self.alnfilename = 'template.aln'

        # create a results file for the Rosetta-generated PDBs to live in
        if not os.path.isdir("results"):
            os.makedirs("results")


    def parseTbl(self):
        # parse tbl.out into infodict
        for line in self.tblhandle:
            if not (line.startswith('#')):

                idandchain = line.split()[0].split(":")
                id = idandchain[0]
                chain = idandchain[1].split(",")[0]

                alnobj = Alignment(id,chain,1,self.fastafile[1].strip(),-1,'',line.split()[8])
                #                                (query)       (template)  score from best domain
                #                                              fill in later
                self.infodict[id] = alnobj
        self.tblhandle.close()

    def parseAln(self):
        # parse aln.out to complete alignment entries in infodict
        for line in self.alnhandle:
            if not (line.startswith('#') or line == '\n'
                    or line == '//\n' or not(line[4]==':')):
                splitline = line.split()
                id = splitline[0].split(":")[0]
                start = splitline[0].split("/")[1].split("-")[0]
                self.infodict[id].tstart = start
                self.infodict[id].talign = splitline[1]
        self.alnhandle.close()




    def writeAlignment(self):
        ## The format string:
        ## Example:
        """
        # score 153.0
        # test 1 LRTYFTTGPQESRAWTINAGMSAPQAAGVIHSDFERGFIRAETIAYKALVEHGSMNAAKEKGLLRSEGKDYVVQEGDVMLFRFN
        # 2OHF_A 305 LEYFFTAGPDEVRAWTIRKGTKAPQAAGKIHTDFEKGFIMAEVMKYEDFKEEGSENAVKAAGKYRQQGRNYIVEDGDIIFFKFN

        # score 152.8
        # test 1 LRTYFTTGPQESRAWTINAGMSAPQAAGVIHSDFERGFIRAETIAYKALVEHGSMNAAKEKGLLRSEGKDYVVQEGDVMLFRFN
        # 2DBY_A 285 -LTFFTAGEKEVRAWTVRRGTKAPRAAGEIHSDMERGFIRAEVIPWDKLVEAGGWARAKERGWVRLEGKDYEVQDGDVIYVLFN

        # score 152.8
        # test 1 LRTYFTTGPQESRAWTINAGMSAPQAAGVIHSDFERGFIRAETIAYKALVEHGSMNAAKEKGLLRSEGKDYVVQEGDVMLFRFN
        # 2DWQ_A 285 -LTFFTAGEKEVRAWTVRRGTKAPRAAGEIHSDMERGFIRAEVIPWDKLVEAGGWARAKERGWVRLEGKDYEVQDGDVIYVLFN

        # score 149.8
        # test 1 LRTYFTTGPQESRAWTINAGMSAPQAAGVIHSDFERGFIRAETIAYKALVEHGSMNAAKEKGLLRSEGKDYVVQEGDVMLFRFN
        # 1JAL_A 279 LQTYFTAGVKEVRAWTVSVGATAPKAAAVIHTDFEKGFIRAEVIAYEDFIQFNGENGAKEAGKWRLEGKDYIVQDGDVMHFRFN

        # score 138.4
        # test 1 LRTYFTTGPQESRAWTINAGMSAPQAAGVIHSDFERGFIRAETIAYKALVEHGSMNAAKEKGLLRSEGKDYVVQEGDVMLFRFN
        # 1NI3_A 308 -INYFTCGEDEVRSWTIRKGTKAPQAAGVIHTDFEKAFVVGEIMHYQDLFDYKTENACRAAGKYLTKGKEYVMESGDIAHWK--
        """

        ## VERY IMPORTANT!!!!!
        ## Remember to uncomment the selected template when you run Rosetta
        ## ALIGNMENT WILL NOT HAPPEN WITHOUT THIS STEP

        __format = string.Template("# score $score\n# $proteinname $qstart $qalign\n# ${pdbid}_${chain} $tstart $talign\n")

        with open(self.alnfilename, 'w') as file:
            for i, key in enumerate(self.infodict.keys()):
                params = {
                    'pdbid': self.infodict[key].pdbid,
                    'chain'  : self.infodict[key].chain,
                    'proteinname': self.protname,
                    'qstart' : self.infodict[key].qstart,
                    'qalign' : self.infodict[key].qalign,
                    'tstart' : self.infodict[key].tstart,
                    'talign' : self.infodict[key].talign,
                    'score' : self.infodict[key].score,
                }
                if (params['tstart'] != -1):

                    file.write('%s\n' % __format.safe_substitute(params))




    def retrievePDB(self):
        # Retrieve template PDBs from RCSB and clean them using clean_pdb.py script
        for key in self.infodict.keys():
            url = "http://www.rcsb.org/pdb/files/"+key+".pdb"
            prefile = key+".pdb"
            if not (os.path.isfile(prefile)):
                wget.download(url,out=prefile)
                clean_pdb.main(prefile,self.infodict[key].chain)






    def createCompModOptions(self,numstructs,db):
        ## Creates the compmod.options file used for minirosetta homology modeling

        begin = ("##INPUT OPTIONS\n"
                "# for minirosetta.linuxgcc\n"
                "\n"
                "# call threading protocol\n"
                "-run:protocol threading\n"
                "\n"
                "-in:file:fullatom\n"
                "-in:file:psipred_ss2 t000_.psipred_ss2\n"
                "-cm:aln_format general\n")

        end = (" ##LOOP BUILDING\n"
                "\n"
                "-loops:frag_sizes 9 3\n"
                "-loops:frag_files aat000_09_05.200_v1_3 aat000_03_05.200_v1_3\n"
                "\n"
                "#use quick CCD to model loops\n"
                "-loops:remodel quick_ccd\n"
                "-idealize_after_loop_close\n"
                "-loops:extended true\n"
                "-loops:build_initial true\n"
                "-loops:relax fastrelax\n"
                "-relax:thorough\n"
                "\n"
                "##OUTPUT OPTIONS\n"
                "-out:overwrite\n"
                "-out:pdb\n"
                "-out:path:pdb results\n"
                "-out:file:fullatom\n")

        compmodstring = begin
        for key in self.infodict.keys():
            compmodstring += "-in:file:template_pdb "+ key + '_' + self.infodict[key].chain +".pdb\n"
        compmodstring += "-in:file:alignment " + self.alnfilename + '\n'
        compmodstring += "-database " + db+ '\n'
        compmodstring += end + '\n'
        compmodstring += "-out:nstruct " + str(numstructs) + '\n'
        compmodstring += "-out:prefix " + self.protname + '\n'

        with open("compmod.options", 'w') as file:
            file.write(compmodstring)



class ManualJM(object):

    def __init__(self,fasta,pdbs,alnpath):
        self.pdbs = pdbs
        self.alnpath = alnpath
        fastaname = fasta.split("/")[-1]
        self.protname = fastaname.split(".")[0]

        if not os.path.isdir("results"):
            os.makedirs("results")

    def retrievePDBs(self):
        with open(self.pdbs,'rb') as pdbs:
            for line in pdbs:
                pdb = line.split()[0]
                chain = line.split()[1]
                url = "http://www.rcsb.org/pdb/files/"+pdb+".pdb"
                prefile = pdb+"pdb"
                file = pdb+chain+".pdb"
                if not (os.path.isfile(file)):
                    wget.download(url,out=prefile)
                    clean_pdb.main(prefile,chain)

    def manual_compModOptions(self,numstructs,db):

        begin = ("##INPUT OPTIONS\n"
                "# for minirosetta.linuxgcc\n"
                "\n"
                "# call threading protocol\n"
                "-run:protocol threading\n"
                "\n"
                "-in:file:fullatom\n"
                "-in:file:psipred_ss2 t000_.psipred_ss2\n"
                "-cm:aln_format general\n")

        end = (" ##LOOP BUILDING\n"
                "\n"
                "-loops:frag_sizes 9 3\n"
                "-loops:frag_files aat000_09_05.200_v1_3 aat000_03_05.200_v1_3\n"
                "\n"
                "#use quick CCD to model loops\n"
                "-loops:remodel quick_ccd\n"
                "-idealize_after_loop_close\n"
                "-loops:extended true\n"
                "-loops:build_initial true\n"
                "-loops:relax fastrelax\n"
                "-relax:thorough\n"
                "\n"
                "##OUTPUT OPTIONS\n"
                "-out:overwrite\n"
                "-out:pdb\n"
                "-out:path:pdb results\n"
                "-out:file:fullatom\n")

        compmodstring = begin
        with open(self.pdbs, 'rb') as pdbs:
            for pdbchain in pdbs:
                p_and_c = pdbchain.split()
                compmodstring += "-in:file:template_pdb "+ p_and_c[0] + '_' + p_and_c[1] +".pdb\n"
        compmodstring += "-in:file:alignment " + self.alnpath + '\n'
        compmodstring += "-database " + db + '\n'
        compmodstring += end + '\n'
        compmodstring += "-out:nstruct " + str(numstructs) + '\n'
        compmodstring += "-out:prefix " + self.protname + '\n'

        with open("compmod.options", 'w') as file:
            file.write(compmodstring)


if __name__ == "__main__":

    print("Running script...")
    t0 = time.clock()

    args = arg_parser()

    if (args.align):
        runjackhmmer(args.protein,args.hmmerpath,args.pdbdb,args.iterations,args.evalue)

    if (args.pdbtemplates and args.alignmentfile):
        jm = ManualJM(args.protein,args.pdbtemplates,args.alignmentfile)
        jm.retrievePDBs()
        jm.manual_compModOptions(args.numdecoys,args.rosettadb)

    elif((args.pdbtemplates and not args.alignmentfile) or (not args.pdbtemplates and args.alignmentfile)):
        print("ERROR: Both PDB templates and alignment file must be given to perform manual setup.")
        sys.exit(1)

    else:
        jm = RunJM(args.protein,'tbl.hits','aln.out')
        jm.parseTbl()
        jm.parseAln()
        jm.writeAlignment()
        jm.retrievePDB()
        jm.createCompModOptions(args.numdecoys,args.rosettadb)

    print("Took", time.clock()-t0, "seconds")