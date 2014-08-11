__author__ = 'kulkarnik'
from Tkinter import *
from lib import annotation
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
import time, sys


def pyplotter2d(finalmat,colors,names):

    # Extract x and y coordinates from finalmat
    x_points = finalmat[:,0]
    y_points = finalmat[:,1]

    # If no colors exist, set default color to blue
    if (colors == []):
        colors = 'b'

    ## Create annotation window with embedded matplotlib
    root = annotation.tk_window_init(x_points,y_points,names,colors)

    # Run Tkinter window
    root.mainloop()


def pyplotter3d(finalmat,colors):

    x_points = finalmat[:,0]
    y_points = finalmat[:,1]
    z_points = finalmat[:,2]

    if (colors == []):
        colors = 'b'

    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.scatter(x_points,y_points,z_points,c=colors)

    plt.show()