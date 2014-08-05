__author__ = 'kulkarnik'
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

def pyplotter2d(finalmat,mdstype,colors):
    if (mdstype == 'm'):
        x_points = finalmat[:,0]
        y_points = finalmat[:,1]
    elif (mdstype == 'c'):
        x_points = []
        y_points = []
        for pair in finalmat:
            end = len(str(pair))
            pairstr = str(pair)[2:end-2].strip().split()
            x_points.append(pairstr[0])
            y_points.append(pairstr[1])

    if (colors == []):
        colors = 'b'

    fig, ax = plt.subplots()
    ax.scatter(x_points, y_points,c=colors)

    plt.show()

def pyplotter3d(finalmat,mdstype,colors):
    if (mdstype == 'm'):
        x_points = finalmat[:,0]
        y_points = finalmat[:,1]
        z_points = finalmat[:,2]
    elif (mdstype == 'c'):
        x_points = []
        y_points = []
        z_points = []
        for triple in finalmat:
            end = len(str(triple))
            tripstr = str(triple)[2:end-2].strip().split()
            x_points.append(tripstr[0])
            y_points.append(tripstr[1])
            z_points.append(tripstr[2])

    if (colors == []):
        colors = 'b'

    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.scatter(x_points,y_points,z_points,c=colors)

    plt.show()