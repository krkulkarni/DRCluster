import approximate_mds
from distance_transform import dist_euclidean
from numpy import array, random
from random import sample

def seedify(dist_mtx,n_seeds):
    seeds = array(sample(dist_mtx,n_seeds))
    return seeds

def perform_nystrom(seed_mtx,dims):
    nystrom_2d = approximate_mds.nystrom(seed_mtx,dims)
