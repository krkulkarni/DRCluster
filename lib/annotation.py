__author__ = 'kulkarnik'

from Tkinter import *
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

class Annotate(object):
    def __init__(self,scatterplot,toolbar,listbox,x_points,y_points,names):
        self.x_points = x_points
        self.y_points = y_points
        self.names = names
        self.selectedpoints = []
        self.scatter = scatterplot
        self.toolbar = toolbar
        self.listbox = listbox

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
        if self.toolbar._active is None:
            self.isPressed = True
            self.x0 = event.xdata
            self.y0 = event.ydata

    def on_release(self, event):

        if self.toolbar._active is None:
            self.isPressed = False
            self.x1 = event.xdata
            self.y1 = event.ydata
            lowerx = min(self.x0,self.x1)
            upperx = max(self.x0,self.x1)
            lowery = min(self.y0,self.y1)
            uppery = max(self.y0,self.y1)

            for i, point in enumerate(self.x_points):
                if (lowerx < self.x_points[i] and self.x_points[i] < upperx and
                            uppery > self.y_points[i] and self.y_points[i] > lowery):
                    if not self.names[i] in self.selectedpoints:
                        self.selectedpoints.append(self.names[i])
                        self.listbox.insert(END,self.names[i])

        self.rect.set_width(0)
        self.rect.set_height(0)
        self.rect.set_xy((self.x0, self.y0))
        self.scatter.figure.canvas.draw()


    def on_pick(self, mouseevent):

        if self.toolbar._active is None:
            self.isPressed = False
            self.ind = mouseevent.ind
            self.x1 = np.take(self.x_points,self.ind)
            self.y1 = np.take(self.y_points,self.ind)
        for i, point in enumerate(self.x_points):
            for x in self.x1:
                for y in self.y1:
                    if (x == self.x_points[i] and y == self.y_points[i]):
                        if not self.names[i] in self.selectedpoints:
                            self.selectedpoints.append(self.names[i])
                            self.listbox.insert(END,self.names[i])


    def on_motion(self,event):
        if self.isPressed:
            self.x1 = event.xdata
            self.y1 = event.ydata
            try:
                self.rect.set_width(self.x1 - self.x0)
                self.rect.set_height(self.y1 - self.y0)
            except TypeError:
                pass
            self.rect.set_xy((self.x0, self.y0))
            self.scatter.figure.canvas.draw()

    def on_key(self,event):
        if event.key == 'escape':
            self.selectedpoints = []
            self.listbox.delete(0, END)

        if event.key == 'r':
            print "Showing all selected points"


def tk_window_init(x_points,y_points,colors,names):
    ## Create Tk main window and assign title and size
    root = Tk()
    root.wm_title("Clustering analysis results")
    root.geometry("1000x700+0+0")

    # Define the quit function used for the Quit Button.
    """# DO NOT USE X TO CLOSE, ONLY USE QUIT"""
    def _quit():
        root.quit()
        root.destroy()

    # Create matplotlib figure and main subplot (ax1)
    fig, ax1 = plt.subplots()
    # Plot points
    ax1.scatter(x_points, y_points,c=colors,picker=2)
    # Remove axes and labels
    ax1.axes.get_xaxis().set_visible(False)
    ax1.axes.get_yaxis().set_visible(False)

    # Create the canvas and place matplotlib figure on canvas
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.show()

    # Create listbox for displaying selected proteins
    listbox = Listbox(master=root,selectmode=EXTENDED)

    # Import toolbar to navigate matplotlib
    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()

    # Create new frame for quit button
    bottomframe = Frame(master=root)

    # Create a scrollbar and attach to listbox
    sb = Scrollbar(master=root,orient=VERTICAL)
    sb.configure(command=listbox.yview)
    listbox.configure(yscrollcommand=sb.set)

    # Pack the components in a reasonable fashion
    toolbar.pack(side=TOP)
    bottomframe.pack(side=BOTTOM)
    canvas.get_tk_widget().pack(side=LEFT, fill=BOTH, expand=1)
    listbox.pack(side=LEFT,fill=BOTH,expand=1,pady=6,padx=3)
    sb.pack(side=LEFT,fill=Y,pady=6,padx=3)

    # Create the point selection and annotation capabilities
    # (More to come soon)
    a = Annotate(ax1,toolbar,listbox,x_points,y_points,names)

    # Create quit button
    quitbutton = Button(master=bottomframe, text='Quit', command=_quit)
    quitbutton.pack(side=LEFT)

    return root
