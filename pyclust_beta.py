import time
import sys
import os
import numpy as np
import json
from lib import updated_readargs, initrun, results_parser, mds_calc, tsne_calc, jsonconv, grouper, plotter

__author__ = "kulkarnik"

##RUN THE CODE
def main():
    ## Get the starting time to measure the time of the run
    print "Running script..."
    t0 = time.clock()



    ## Obtain all the info from the arguments passed
    args = updated_readargs.arg_parser()
    print "Parsed arguments"


    ## Define the paths of BLAST results file and names file
    resultsfile = args.directory + "/results.out"
    fastafile = args.directory + "/" + args.directory.split('/')[-1] + ".fas"
    matrixpath = args.directory + "/temp/mds.hdf5"
    coordspath = args.directory + "/temp/" + args.directory.split('/')[-1] + "_"+ args.type + "_" + "coords.npy"
    jsonpath = args.directory + "/temp/file.json"
    inity = args.directory + "/temp/inity.npy"
    grouppath = args.directory +"/temp/groups.txt"

    if (args.reinitialize):
        try:
            os.remove(inity)
        except OSError:
            pass

    ## Obtain colors and names from names file
    names,colors,lines,seqs = initrun.read_fasta(fastafile,args.color)
    print "Read FASTA file"
    ## Obtain the handle to the results file
    tabParser, tabHandle = initrun.open_file(resultsfile)
    print "Opened file"

    ## Check if coordinates have already been mapped
    if (args.load):
        print "Loading coordinates from path"
        matrix = np.load(args.load)

    elif (args.parse):
        hdfmat, mathandle = initrun.create_matrix(args.value,names,matrixpath)
        print "Initialized matrix"

        print "Parsing results"
        results_parser.next_line_original_format(args.value,tabParser,tabHandle,names,hdfmat)

    else:
        print "Retrieving matrix"
        hdfmat, mathandle = initrun.get_matrix(matrixpath)

        ## Run the appropriate dimensionality reduction algorithm
        ## -mdsonly = metric MDS with sklearn's manifold package
        ## -snemds = preprocess to "points/10" dimensions with MDS, then t-SNE for reduced matrix
        ## -snepca = preprocess to "points/10" dimensions with PCA, then t-SNE for reduced matrix

    if (args.cluster):
        if (args.type == "mdsonly"):
            print "Performing MDS"
            matrix = mds_calc.metric_mds(hdfmat,int(args.dimension))

        elif (args.type == "snemds"):
            print "Performing t-SNE with MDS preprocessing"

            ## Partially reduce dimensionality of HDF5 matrix to 1/10th of original size or maximum of 400
            tempred = min(int(len(names)/10),400)
            print "Preprocessing the data using MDS..."
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
        mathandle.close()

    else:
        print "Loading coordinates from temp: " + args.type
        matrix = np.load(coordspath)
        mathandle.close()

    if (args.group):
        #groupids = grouping.findkmeans(matrix)
        groupids = grouper.findgroups(matrix,args.group)
        np.savetxt(grouppath,groupids)
        if (args.color == "group"):
            colors = groupids
    elif (args.color == 'group'):
        groupids = np.loadtxt(grouppath)
        if (args.color == 'group'):
            colors = groupids

    #convert colors to gradient colors
    #colors = coloring.convert(colors)

    print "Took", time.clock()-t0, "seconds"

    #with open(jsonpath, 'w') as jsonout:
    #    json.dump(jsonconv.jsonmaker(colors,lines,matrix,args.format), jsonout, indent=2)

    ## Plot the results with matplotlib's PyPlot
    if (args.plot):
        print "Plotting",len(names), "points"
        if (int(args.dimension) == 2):
            plotter.pyplotter2d(matrix,colors,names,seqs,args.directory)

        elif (int(args.dimension) == 3):
            plotter.pyplotter3d(matrix,colors)

if (__name__ == '__main__'):
    main()






