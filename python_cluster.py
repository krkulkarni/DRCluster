import time, sys, os
import numpy as np
from Tkinter import *
from lib import readargs, initrun, results_parser, mds_calc, tsne_calc, plotter



__author__ = "kulkarnik"


##RUN THE CODE
if (__name__ == '__main__'):

    ## Get the starting time to measure the time of the run
    print "Running script..."
    t0 = time.clock()



    ## Obtain all the info from the arguments passed
    args = readargs.arg_parser()
    print "Parsed arguments"


    """
    If necessary, run BLAST locally (try to run on a faster machine instead)
    #initrun.run_blast_tab(queryname,dbname,outfile,fmt,dbsize,ecutoff)
    """


    ## Define the paths of BLAST results file and names file
    resultsfile = args.directory + "/results.out"
    fastafile = args.directory + "/names_all"
    matrixpath = args.directory + "/temp/mds.hdf5"
    coordspath = args.directory + "/temp/coords.npy"
    inity = args.directory + "/temp/inity.npy"

    if (args.reinitialize):
        try:
            os.remove(inity)
        except OSError:
            pass

    ## Obtain colors and names from names file
    names,colors = initrun.read_fasta(fastafile,args.format)
    print "Read FASTA file"
    ## Obtain the handle to the results file
    tabParser, tabHandle = initrun.open_file(resultsfile)
    print "Opened file"

    if args.precoordinated == True:
        print "Loading coordinates"
        matrix = np.load(coordspath)
    else:
        ## Check if results are preparsed
        ##
        ## If no, create an HDF5 formatted matrix and initialize
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



        ## Run the appropriate dimensionality reduction algorithm
        ## -mdsonly = metric MDS with sklearn's manifold package
        ## -snemds = preprocess to "points/10" dimensions with MDS, then t-SNE for reduced matrix
        ## -snepca = preprocess to "points/10" dimensions with PCA, then t-SNE for reduced matrix

        ## (still working on -n = nystrom MDS with pycogent's approximate_mds package)
        if (args.type == "mdsonly"):
            print "Performing MDS"
            matrix = mds_calc.metric_mds(hdfmat,int(args.dimension))

        elif (args.type == "snemds"):
            print "Performing t-SNE with MDS preprocessing"

            ## Partially reduce dimensionality of HDF5 matrix to 1/10th of original size or maximum of 400
            tempred = min(int(len(names)/10),400)
            print "Reducing to", tempred, "dimensions"

            tempmatrix = mds_calc.metric_mds(hdfmat,tempred)
            matrix = tsne_calc.tsne(inity,False,tempmatrix,no_dims=int(args.dimension),initial_dims=tempred)

        elif (args.type == "snepca"):
            print "Performing t-SNE with PCA preprocessing"

            ## Partially reduce dimensionality of HDF5 matrix to 1/10th of original size or maximum of 400
            tempred = min(int(len(names)/10),400)
            print "Reducing to", tempred, "dimensions"

            matrix = tsne_calc.tsne(inity,True,hdfmat,no_dims=int(args.dimension),initial_dims=tempred)

        elif (args.type == "sneonly"):
            print "CAUTION: Performing t-SNE on full dissimilarity matrix"

            if (len(names) > 2000):
                print "Too many proteins to perform t-SNE directly"
                sys.exit(2)
            matrix = tsne_calc.tsne(inity,False,hdfmat[...],no_dims=int(args.dimension),initial_dims=len(names))

        # model = TSNE(n_components=3,metric="precomputed")
        # matrix = model.fit_transform(hdfmat)
        #else: matrix = mds_calc.nystrom_frontend(len(names),math.sqrt(len(names)),2,mds_calc.getdist,hdfmat)

    # save coordinates to file
    np.save(coordspath,matrix)
    print "Took", time.clock()-t0, "seconds"


    ## Plot the results with matplotlib's PyPlot
    if (args.plot):
        print "Plotting",len(names), "points"

        if (int(args.dimension) == 2):
            plotter.pyplotter2d(matrix,colors,names)

        elif (int(args.dimension) == 3):
            plotter.pyplotter3d(matrix,colors)






