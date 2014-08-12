__author__ = 'kulkarnik'
from Tkinter import *
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from lib import annotation

from mpl_toolkits.mplot3d import Axes3D
import time, sys


def pyplotter2d(finalmat,colors,names):

    # Extract x and y coordinates from finalmat
    x_points = finalmat[:,0]
    y_points = finalmat[:,1]

    # If no colors exist, set default color to blue
    if (colors == []):
        colors = 'b'

    # Create Tk window and handles to subplot, toolbar and listbox
    root,ax1,toolbar,listbox = annotation.tk_window_init(x_points,y_points,names,colors)

    # Create capability to select points
    a = annotation.Annotate(ax1,toolbar,listbox,x_points,y_points,names,colors)

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