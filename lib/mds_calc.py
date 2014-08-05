__author__ = 'kulkarnik'
from lib import approximate_mds
from sklearn import manifold
import numpy as np

##perform the MDS on HDF5 matrix
def cmds(matrix,dim,coordsfile):
    if (dim==2):
        ## MDS with 2 dimensions (default)
        (newmat,eig) = approximate_mds.cmds_tzeng(matrix,dim=2)

        np.save(coordsfile,newmat)
        return newmat
    elif (dim==3):
        ## MDS with 3 dimensions
        (newmat,eig) = approximate_mds.cmds_tzeng(matrix,dim=3)
        ## call plotter with 3D graph
        return newmat

def metric_mds(mat,dim,coordsfile):
    mds = manifold.MDS(n_components=dim, metric=True,
                           max_iter=300,dissimilarity="precomputed",
                           n_jobs=-1,n_init=1,random_state=1)
    coords = mds.fit(mat).embedding_

    return coords

# def nystrom_frontend(num_objects, num_seeds, dim, dist_func,permute_order=True):
#     (seed_mat,restore_ids) = approximate_mds.build_seed_matrix(num_objects,num_seeds,dist_func,permute_order)
#     mdscoords = approximate_mds.nystrom(seed_mat,dim)
#     return mdscoords[restore_ids]
#
# def getdist(i,j,hdfmat):
#     return hdfmat[i,j]
