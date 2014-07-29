from lib import readargs, initrun, results_parser, mds_calc, plotter
import time

##
##GLOBAL VARIABLES
##

__author__ = "kulkarnik"

resultsfile = "/Users/kulkarnik/Research/MDSCluster_2014/Ub_diff/curr/results.out"
fastafile = "/Users/kulkarnik/Research/MDSCluster_2014/Ub_diff/curr/names_ub"


##RUN THE CODE
if (__name__ == '__main__'):
    t0 = time.clock()
    print "Running script..."
    bitOrE, dim, mdstype, fmat = readargs.arg_parser()
    print "Parsed arguments"

    #run_blast_tab(queryname,dbname,outfile,fmt,dbsize,ecutoff)
    names,colors = initrun.read_fasta(fastafile)
    print "Read FASTA file"

    tabParser, tabHandle = initrun.open_file(resultsfile)
    print "Opened file"

    hdfmat = initrun.create_matrix(bitOrE,names)
    print "Initialized matrix"

    print "Parsing results"
    if (fmat == 'mod'):
        results_parser.next_line_modified_format(bitOrE,tabParser,tabHandle,names,hdfmat)
    else:
        results_parser.next_line_original_format(bitOrE,tabParser,tabHandle,names,hdfmat)


    print "Performing MDS"
    if (mdstype == "m"):
        matrix = mds_calc.metric_mds(hdfmat)
    elif (mdstype == "c"):
        matrix = mds_calc.cmds(hdfmat,dim)

    #matrix = mds_calc.nystrom_frontend(len(names),math.sqrt(len(names)),2,mds_calc.getdist)

    plotter.pyplotter2d(matrix,mdstype,colors,names,t0)





