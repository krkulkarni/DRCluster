__author__ = 'kulkarnik'
import numpy as np

# inityfull = np.load("/Users/kulkarnik/Research/MDSCluster_2014/pyclust/full_ub/temp/inity.npy")
# np.savetxt("/Users/kulkarnik/Research/MDSCluster_2014/pyclust/full_ub/temp/inity.txt",inityfull)
#
# initysub = np.load("/Users/kulkarnik/Research/MDSCluster_2014/pyclust/reps/temp/inity.npy")
# np.savetxt("/Users/kulkarnik/Research/MDSCluster_2014/pyclust/reps/temp/inity.txt",initysub)

inityfinal = np.loadtxt("/Users/kulkarnik/Research/MDSCluster_2014/pyclust/full_ub/temp/inity.txt")
np.save("/Users/kulkarnik/Research/MDSCluster_2014/pyclust/total_ub/temp/inity.npy",inityfinal)