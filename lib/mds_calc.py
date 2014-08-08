__author__ = 'kulkarnik'

from sklearn import manifold

## Perform the MDS on HDF5 matrix
## n_components = final number of dimensions
## n_jobs = number of cores to use in calculation (-1 means use all cores)

def metric_mds(mat,dim):
    print "Preprocessing the data using MDS..."
    mds = manifold.MDS(n_components=dim, metric=True,
                           max_iter=300,dissimilarity="precomputed",
                           n_jobs=1,n_init=1,random_state=1)
    coords = mds.fit(mat).embedding_

    return coords

# def nystrom_frontend(num_objects, num_seeds, dim, dist_func,permute_order=True):
#     (seed_mat,restore_ids) = approximate_mds.build_seed_matrix(num_objects,num_seeds,dist_func,permute_order)
#     mdscoords = approximate_mds.nystrom(seed_mat,dim)
#     return mdscoords[restore_ids]
#
# def getdist(i,j,hdfmat):
#     return hdfmat[i,j]
