import os, sys
import numpy as np
import tsne
from lib import mds_calc, tsne_calc

__author__ = 'kulkarnik'
## The Algorithm class holds the three methods of dimensionality reduction.

## SVD - tSNE hybrid: svdsne()
## First performs singular value decomposition directly on sparse matrix
## to reduce the noise of dataset.
## Then, runs t-SNE with Barnes-Hut approximation on intermediate dataset
## to reduce to coordinate embedding

## Euclidean Multidimensional Scaling: mdsonly()
## Use with caution! It fails on large datasets
## Convert sparse matrix to dense matrix
## and perform linear dimensionality reduction
## treating similarity values as Euclidean distances

## t-SNE only: sneonly()
## Use with caution! This method fails under most conditions,
## due to the organization of the sparse similarity matrix
## Runs t-SNE algorithm directly on full matrix, without intermediate reduction,
## and without Barnes-Hut approximation

class Algorithm(object):

    def __init__(self,scipymat,points,dim,roadmap):
        self.scipymat = scipymat
        self.points = points
        self.pointslen = len(points)
        self.dim = dim
        dirname = os.path.dirname(os.path.realpath(__file__))
        self.seed = self._generate_seed_matrix()
        # numpyfile = "{}/seed.txt".format(dirname)
        # if roadmap:
        #     self.seed = self._normalize_seed(np.loadtxt(roadmap),numpyfile)
        # else:
        #     self.seed = self._normalize_seed(np.genfromtxt(numpyfile,skip_footer=120000-self.pointslen))

    def coordsgen(self,points):
        amino_acids = {
            "G": -.7444439610016066,
            "A": -.8120238054672788,
            "V": -.6814501991655911,
            "L": -.4655948851103902,
            "I": -.522471902925321,
            "P": -.9733493722846006,

            "F": -.12813674753375167,
            "Y": -.2831025840191209,
            "W": -.2060643234971105,

            "S": .22299402164794557,
            "T": .3971822974996705,
            "C": .06914052935419317,
            "N": .33157208503360835,
            "M": .4315294239386928,
            "Q": .19041624661608026,

            "K": .6800921873765579,
            "R": .7001202297117273,
            "H": .7478338718290908,

            "D": .8803319774576138,
            "E": .9738805050045272
        }

        for key, point in points.iteritems():
            yield sum(amino_acids[aa] for aa in point.seq[0:10]), \
                  sum(amino_acids[aa] for aa in point.seq[-10:-1])

    def _generate_seed_matrix(self):
        seed = []
        for x,y in self.coordsgen(self.points):
            seed.append([x,y])
        return np.array(seed)

    # def _normalize_seed(self,roadmap,*args):
    #     mappoints, _ = roadmap.shape
    #     if (mappoints == self.pointslen):
    #         print("Roadmap returned")
    #         return roadmap
    #
    #     elif (mappoints < self.pointslen):
    #         print("Added points to roadmap")
    #         newpoints = self.pointslen - mappoints
    #         newseed = np.genfromtxt(args[0],skip_header=mappoints,skip_footer=120000-self.pointslen)
    #         #newzeroseed = np.zeros((newpoints,self.dim),dtype=np.float64)
    #         try:
    #             return np.concatenate((roadmap,newseed))
    #         except ValueError:
    #             # Only happens if newseed has one member
    #             return np.concatenate((roadmap,[newseed]))
    #     else:
    #         # If mappoints > self.pointslen
    #         return roadmap[:self.pointslen]

    def svdsne(self,perp,theta):
        print("Performing svdsne")
        tempred = min(self.pointslen/10,50)
        print("Reducing to {} dimensions with SVD".format(tempred))
        tempmatrix = mds_calc.svd(self.scipymat,tempred)
        matrix = tsne.bh_sne(tempmatrix,self.seed,perplexity=perp,theta=theta)
        return matrix


    def mdsonly(self,*args):
        print("Performing mdsonly")
        if (self.pointslen > 2000):
            print "Too many proteins to perform MDS directly"
            sys.exit(2)
        matrix = mds_calc.metric_mds(self.scipymat.toarray(),self.dim)
        return matrix


    def sneonly(self,reinit,directory,**kwargs):
        print("Performing sneonly")
        inity = "{}/inity.npy".format(directory)
        if (reinit):
            try:
                os.remove(inity)
            except OSError:
                print "No initial inity file"
        if (self.pointslen > 2000):
            print "Too many proteins to perform t-SNE directly"
            sys.exit(2)
        matrix = tsne_calc.tsne(inity,False,self.scipymat.toarray(),
                                no_dims=self.dim,
                                initial_dims=self.pointslen)
        return matrix
