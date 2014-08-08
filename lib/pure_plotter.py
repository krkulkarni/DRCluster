__author__ = 'kulkarnik'
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from Tkinter import *

from mpl_toolkits.mplot3d import Axes3D
import time, sys



class Annotate(object):
    def __init__(self,scatterplot,top,x_points,y_points,names):
        self.x_points = x_points
        self.y_points = y_points
        self.names = names
        self.selectedpoints = []
        self.scatter = scatterplot
        self.top = top

        self.lb1 = Listbox(self.top)
        self.lb1.pack()
        # try:
        #     self.top.mainloop()
        # except:
        #     sys.exit(0)

        # self.data = []
        # for i in xrange(27):
        #     self.data.append([""])
        # self.table = self.tableplot.table(cellText=self.data,loc='center')

        self.rect = Rectangle((0,0), 0, 0,facecolor='grey', alpha=0.3)
        self.x0 = 0
        self.y0 = 0
        self.x1 = 0
        self.y1 = 0
        self.isPressed = False
        self.scatter.add_patch(self.rect)

        self.scatter.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.scatter.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.scatter.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.scatter.figure.canvas.mpl_connect('pick_event', self.on_pick)
        self.scatter.figure.canvas.mpl_connect('key_press_event', self.on_key)

    def on_press(self, event):
        if self.scatter.figure.canvas.manager.toolbar._active is None:
            self.isPressed = True
            self.x0 = event.xdata
            self.y0 = event.ydata

    def on_release(self, event):

        if self.scatter.figure.canvas.manager.toolbar._active is None:
            self.isPressed = False
            self.x1 = event.xdata
            self.y1 = event.ydata

        for i, point in enumerate(self.x_points):
            if (self.x0 < self.x_points[i] and self.x_points[i] < self.x1 and
                        self.y0 > self.y_points[i] and self.y_points[i] > self.y1):
                if not self.names[i] in self.selectedpoints:
                    print self.names[i]
                    self.selectedpoints.append(self.names[i])
                    self.lb1.insert(END,self.names[i])

        self.rect.set_width(0)
        self.rect.set_height(0)
        self.rect.set_xy((self.x0, self.y0))
        self.scatter.figure.canvas.draw()


    def on_pick(self, mouseevent):

        if self.scatter.figure.canvas.manager.toolbar._active is None:
            self.isPressed = False
            self.ind = mouseevent.ind
            self.x1 = np.take(self.x_points,self.ind)
            self.y1 = np.take(self.y_points,self.ind)
        for i, point in enumerate(self.x_points):
            for x in self.x1:
                for y in self.y1:
                    if (x == self.x_points[i] and y == self.y_points[i]):
                        if not self.names[i] in self.selectedpoints:
                            print self.names[i]
                            self.selectedpoints.append(self.names[i])
                            self.lb1.insert(END,self.names[i])


    def on_motion(self,event):
        if self.isPressed:
            self.x1 = event.xdata
            self.y1 = event.ydata
            self.rect.set_width(self.x1 - self.x0)
            self.rect.set_height(self.y1 - self.y0)
            self.rect.set_xy((self.x0, self.y0))
            self.scatter.figure.canvas.draw()

    def on_key(self,event):
        if event.key == 'escape':
            self.selectedpoints = []
            print "\nCleared selected points"

        if event.key == 'r':
            print "\nShowing all selected points"


def pyplotter2d(finalmat,colors,names,root):

    x_points = finalmat[:,0]
    y_points = finalmat[:,1]

    if (colors == []):
        colors = 'b'

    fig, ax1 = plt.subplots()
    ax1.scatter(x_points, y_points,c=colors,picker=2)


    a = Annotate(ax1,root,x_points,y_points,names)


    plt.show()
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