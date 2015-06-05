#!/usr/bin/env python

import time
import sys
import os
import numpy as np
from lib import updated_readargs, initrun, results_parser, mds_calc, tsne_calc, grouper, plotter
from scipy import sparse
import tsne
#import json
#from lib import jsonconv

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
    hmmerfile = args.directory + "/tbl.hits"
    fastafile = args.directory + "/" + args.directory.split('/')[-1] + ".fas"
    matrixpath = args.directory + "/temp/data.txt"
    coordspath = args.directory + "/temp/" + args.directory.split('/')[-1] + "_"+ args.type + "_" + "coords.npy"
    jsonpath = args.directory + "/temp/file.json"
    inity = args.directory + "/temp/inity.npy"
    grouppath = args.directory +"/temp/groups.txt"
    sizepath = args.directory + "/temp/sizes.txt"
    reppath = args.directory +"/temp/reps.txt"

    if (args.reinitialize):
        try:
            os.remove(inity)
        except OSError:
            print "No initial inity file"

    ## Obtain colors and names from names file
    points = initrun.read_fasta(fastafile)
    print "Read FASTA file"

    ## Check if coordinates have already been mapped
    if (args.load):
        print "Loading coordinates from path"
        matrix = np.load(args.load)

    elif (args.parse):

        ## Obtain the handle to the results file
        if (args.search == 'blast'):
            tabParser, tabHandle = initrun.open_file(resultsfile)
        elif (args.search == "hmmer"):
            tabParser, tabHandle = initrun.open_file(hmmerfile)
        print "Opened BLAST results file"

        # hdfmat, mathandle = initrun.create_matrix(args.value,points,matrixpath)
        # print "Initialized matrix"

        print "Parsing results"
        row,col,data = results_parser.next_line_original_format(args.value,tabParser,tabHandle,points,args.search)
        savemat = np.vstack((row,col,data))
        np.savetxt(matrixpath,savemat)

        scipymat = sparse.coo_matrix((data,(row,col)),shape=(len(points),len(points)))

    else:
        print "Retrieving matrix"
        savemat = initrun.get_matrix(matrixpath)
        row = savemat[0]
        col = savemat[1]
        data = savemat[2]

        scipymat = sparse.coo_matrix((data,(row,col)),shape=(len(points),len(points)))
        ## Run the appropriate dimensionality reduction algorithm
        ## -mdsonly = metric MDS with sklearn's manifold package
        ## -svdsne = preprocess to "points/10" dimensions with MDS, then t-SNE for reduced matrix
        ## -sneonly = t-SNE with van der Maatan algorithm

    if (args.cluster):
        if (args.type == "mdsonly"):
            print "Performing MDS"
            if (len(points) > 2000):
                print "Too many proteins to perform MDS directly"
                sys.exit(2)
            matrix = mds_calc.metric_mds(scipymat.toarray(),int(args.dimension))

        elif (args.type == "svdsne"):
            print "Performing t-SNE with MDS preprocessing"

            ## Partially reduce dimensionality of HDF5 matrix to 1/10th of original size or maximum of 400
            tempred = min(int(len(points)/10),400)
            print "Preprocessing the data using SVD..."
            print "Reducing to", tempred, "dimensions"

            tempmatrix = mds_calc.svd(scipymat,tempred)
            matrix = tsne.bh_sne(tempmatrix)
            # matrix = tsne_calc.tsne(inity,False,tempmatrix,no_dims=int(args.dimension),initial_dims=tempred)

        elif (args.type == "sneonly"):
            print "CAUTION: Performing t-SNE on full dissimilarity matrix"

            if (len(points) > 2000):
                print "Too many proteins to perform t-SNE directly"
                sys.exit(2)
            matrix = tsne_calc.tsne(inity,False,scipymat.toarray(),no_dims=int(args.dimension),initial_dims=len(points))
        else:
            print "Some error occurred! Please try again."
            sys.exit(2)
            # model = TSNE(n_components=3,metric="precomputed")
            # matrix = model.fit_transform(hdfmat)
            #else: matrix = mds_calc.nystrom_frontend(len(names),math.sqrt(len(names)),2,mds_calc.getdist,hdfmat)

        # save coordinates to file
        np.save(coordspath,matrix)

    else:
        try:
            print "Loading coordinates from temp: " + args.type
            matrix = np.load(coordspath)
        except IOError:
            print "No coordinates in temp"

    if (args.group):
        labels = grouper.dbscan(matrix)
        #groupids = grouping.findkmeans(matrix)
        # print "Grouping points"
        # groupids, reps, sizes = grouper.findgroups(matrix,args.group,points)
        # np.savetxt(grouppath,groupids)
        # np.savetxt(sizepath, sizes)
        # grouper.savereps(reppath,reps,args.group)

    # select correct color array
    colors = np.zeros(shape=(len(points)))
    if (args.color == 'pfam'):
        for point in points:
            colors[points[point].index] = points[point].pfamnum
    elif (args.color == 'mod'):
        for point in points:
            colors[points[point].index] = points[point].modcolor
    elif (args.color == 'group'):
        colors = labels
        # colors = np.loadtxt(grouppath)
    else:
        print "Some error occurred"
        sys.exit(1)

    print len(matrix), "points in dataset"
    print "Took", time.clock()-t0, "seconds"

    #with open(jsonpath, 'w') as jsonout:
    #    json.dump(jsonconv.jsonmaker(colors,lines,matrix,args.format), jsonout, indent=2)

    ## Plot the results with matplotlib's PyPlot
    if (args.plot):
        print "Plotting",len(matrix), "points"
        ## load sizes manually!!!
        if (args.directory.split('/')[-1] == "total_ub" or args.directory.split('/')[-1] == "full_ub"):
            point_sizes = np.loadtxt(args.directory + "/allsizes.txt")
        else:
            point_sizes = 20
        if (int(args.dimension) == 2):
            plotter.pyplotter2d(matrix,colors,args.directory,points,point_sizes)

        elif (int(args.dimension) == 3):
            plotter.pyplotter3d(matrix,colors)

if (__name__ == '__main__'):
    main()






