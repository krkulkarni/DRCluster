import os
import sys
sys.path.append("../../lib/tsne/")
import numpy as np
import tsne

import mds_calc, tsne_calc

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

    def __init__(self,scipymat,points,dim, seedpath):
        self.scipymat = scipymat
        self.points = points
        self.pointslen = len(points)
        self.dim = dim
        if (seedpath):
            self.seed = np.loadtxt(seedpath)
        else:
            self.seed = self._generate_seed_matrix()


    def _coordsgen(self,points):
        amino_acids = {
            "G": -.7444439610016066, #
            "A": -.8120238054672788, # ALL
            "V": -.6814501991655911, # NONPOLAR
            "L": -.4655948851103902, # RESIDUES
            "I": -.522471902925321,  # (-1.0,-0.4)
            "P": -.9733493722846006, #

            "F": -.12813674753375167,# ALL AROMATIC
            "Y": -.2831025840191209, # RESIDUES
            "W": -.2060643234971105, # (-0.4,-0.1)

            "S": .22299402164794557, #
            "T": .3971822974996705,  # ALL POLAR,
            "C": .06914052935419317, # UNCHARGED
            "N": .33157208503360835, # RESIDUES
            "M": .4315294239386928,  # (-0.1,0.5)
            "Q": .19041624661608026, #

            "K": .6800921873765579,  # ALL POLAR,
            "R": .7001202297117273,  # POSITIVE RESIDUES
            "H": .7478338718290908,  # (0.5,0.8)

            "D": .8803319774576138,  # ALL POLAR, NEGATIVE
            "E": .9738805050045272,   # RESIDUES
                                     # (0.8,1.0)

            "X": 0.0                 # Random residue
        }

        for i, (key, point) in enumerate(points.iteritems(),start=1):
            if (len(point.seq) >= 10):
                x_coord = sum(amino_acids[aa] for aa in point.seq[0:10])
                y_coord = sum(amino_acids[aa] for aa in point.seq[-10:-1])
            else:
                x_coord = sum(amino_acids[aa] for aa in point.seq[0:len(point.seq)-1])
                y_coord = sum(amino_acids[aa] for aa in point.seq[0:len(point.seq)-1])

            if (len(point.seq) >= 20):
                z_coord = sum(amino_acids[aa] for aa in point.seq[0:20:2])
            else:
                z_coord = sum(amino_acids[aa] for aa in point.seq[0:len(point.seq)-1:2])
            yield i,\
                  x_coord, y_coord, z_coord

    def _generate_seed_matrix(self):
        seed = []
        print("Calculating initial seeding matrix")
        print("Calculating coordinates for point {} of {}".format("0",self.pointslen))
        for i,x,y,z in self._coordsgen(self.points):
            if (self.dim == 2):
                seed.append([x,y])
            elif (self.dim == 3):
                seed.append([x,y,z])
            elif (self.dim == 4):
                seed.append([x,y,z,0])
            else:
                print("Invalid dimensionality! Strange error")
                raise NotImplementedError
            if (i%1000 == 0):
                print("Calculating coordinates for point {} of {}".format(i,self.pointslen))
        return np.array(seed)


    def svdsne(self,perp,theta,maxiter):
        print("Performing svdsne")
        tempred = min(self.pointslen/5,200)
        print("Reducing to {} dimensions with SVD".format(tempred))
        tempmatrix = mds_calc.svd(self.scipymat,tempred)
        matrix = tsne.bh_sne(tempmatrix,self.seed,maxiter,
                             d=self.dim,perplexity=perp,
                             theta=theta,pca_d=None)
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
