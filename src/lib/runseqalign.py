__author__ = 'kulkarnik'

import subprocess

class Align(object):

    def __init__(self,fastapath,exe,directory):

        self.path = fastapath
        self.base = fastapath.split("/")[-1].split(".")[0]
        self.exe = exe
        self.directory = directory


    def _execute_commands(self,command_array, wait=True):
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

    def runjackhmmer(self,iter,evalue):

        ## This function accepts the path to the target fasta sequence and the fasta PDB database
        ## and creates a domain hits file (dom.hits),
        ## table of hits file (tbl.hits)
        ## jackhmmer alignment file (aln.out)
        ## and a hmmer log file (hmmerlog)
        ## see Documentation at http://hmmer.janelia.org/software for additional info

        tblout = "{}/{}.jackhmmer_tbl".format(self.directory,self.base)

        comm = [
             [ self.exe, '--tblout', tblout,
               '-E', str(evalue), '-N', str(iter),
               self.path, self.path,
             ]
           ]

        print(' '.join(comm[0]))
        return self._execute_commands(comm)
        ## Returns jackhmmer ouput

    def runblast(self):

        results = "{}/{}.blast_tbl".format(self.directory)
        db = "{}/{}".format(self.directory,self.base)

        comm = [
             [ self.exe, '-query', self.path,
               '-db', db, '-outfmt', 6, '-out', results
             ]
           ]

        print(' '.join(comm[0]))
        r = self._execute_commands(comm)
        #Returns Blast output (should be nothing)