from lib import readargs, initrun, results_parser, mds_calc, plotter
import time, sys
import numpy as np
import tsne_calc
from sklearn.manifold import TSNE



##
##GLOBAL VARIABLES
##

__author__ = "kulkarnik"

## Define the paths of BLAST results file and names file


##RUN THE CODE
if (__name__ == '__main__'):

    ## Get the starting time to measure the time of the run
    print "Running script..."
    t0 = time.clock()



    ## Obtain all the info from the arguments passed
    args = readargs.arg_parser()
    print "Parsed arguments"


    resultsfile = "/Users/kulkarnik/Research/MDSCluster_2014/BLAST+/"+args.directory+"/results.out"
    fastafile = "/Users/kulkarnik/Research/MDSCluster_2014/BLAST+/"+args.directory+"/names_all"
    matrixpath = "/Users/kulkarnik/Research/MDSCluster_2014/BLAST+/"+args.directory+"/temp/mds.hdf5"
    coordspath = "/Users/kulkarnik/Research/MDSCluster_2014/BLAST+/"+args.directory+"/temp/coords.npy"
    inity = "/Users/kulkarnik/Research/MDSCluster_2014/BLAST+/"+args.directory+"/temp/inity.npy"

    """
    If necessary, run BLAST locally (try to run on a faster machine instead)
    #run_blast_tab(queryname,dbname,outfile,fmt,dbsize,ecutoff)
    """



    ## Obtain colors and names from names file
    names,colors = initrun.read_fasta(fastafile)
    print "Read FASTA file"
    ## Obtain the handle to the results file
    tabParser, tabHandle = initrun.open_file(resultsfile)
    print "Opened file"



    ## Check if results are preparsed
    ##
    ## If not, create an HDF5 formatted matrix and initialize
    ## Then, parse BLAST results and populate the HDF5 matrix
    ##
    ## If yes, then obtain populated matrix from file

    if args.preparsed:
        print "Retrieving matrix"
        hdfmat = initrun.get_matrix(matrixpath)
    else:
        hdfmat = initrun.create_matrix(args.value,names,matrixpath)
        print "Initialized matrix"

        print "Parsing results"
        if (args.format == 'mod'):
            results_parser.next_line_modified_format(args.value,tabParser,tabHandle,names,hdfmat)
        elif (args.format == 'orig'):
            results_parser.next_line_original_format(args.value,tabParser,tabHandle,names,hdfmat)

    ## Check if coordinates have already been mapped
    if args.precoordinated == True:
        print "Loading coordinates"
        matrix = np.load(coordspath)

    else:
        ## Run the appropriate MDS algorithm
        ## -m = metric MDS with sklearn's manifold package
        ## -c = classical MDS with pycogent's approximate_mds package
        ## (still working on -n = nystrom MDS with pycogent's approximate_mds package)
        print "Performing MDS"
        if (args.mdstype == "m"):
            tempmatrix = mds_calc.metric_mds(hdfmat,50,coordspath)

        elif (args.mdstype == "c"):
            matrix = mds_calc.cmds(hdfmat,args.dimension,coordspath)

        print "Performing t-SNE"
        matrix = tsne_calc.tsne(inity,tempmatrix,no_dims=3,initial_dims=50)
        # model = TSNE(n_components=3,metric="precomputed")
        # matrix = model.fit_transform(hdfmat)
        #else: matrix = mds_calc.nystrom_frontend(len(names),math.sqrt(len(names)),2,mds_calc.getdist,hdfmat)

    np.save(coordspath,matrix)
    print "Took", time.clock()-t0, "seconds"


    ## Plot the results with matplotlib's pyplot
    if (args.plot):
        print "Plotting",len(names), "points"
        if (int(args.dimension) == 2):
            plotter.pyplotter2d(matrix,args.mdstype,colors)

        elif (int(args.dimension) == 3):
            plotter.pyplotter3d(matrix,args.mdstype,colors)





