__author__ = 'kulkarnik'

from Tkinter import *
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Rectangle
import numpy as np
import csv

class Annotate(object):
    def __init__(self,scatterplot,toolbar,listbox,infobox,x_points,y_points,colors,dir,points,sizes):
        self.path = dir + '/selected.fas'
        self.x_points = x_points
        self.y_points = y_points
        self.selectedpoints = []
        self.selectedseqs = []
        self.scatter = scatterplot
        self.toolbar = toolbar
        self.listbox = listbox
        self.infobox = infobox
        self.colors = colors
        self.points = points
        self.sizes = sizes
        self.cmaps = [cm.prism,cm.Accent,cm.gist_ncar,cm.Paired,cm.rainbow,cm.Blues]
        self.cmapindex = 0

        self.scatter.scatter(self.x_points, self.y_points,c=self.colors,cmap=self.cmaps[self.cmapindex],s=self.sizes)

        self.rect = Rectangle((0,0), 0, 0,facecolor='grey', alpha=0.3)
        self.scatter.add_patch(self.rect)
        self.x0 = 0
        self.y0 = 0
        self.x1 = 0
        self.y1 = 0
        self.isPressed = False

        self.scatter.figure.canvas.draw()
        self.connect()

    def onselect(self, event):
    # Note here that Tkinter passes an event object to onselect()
        w = event.widget
        index = (w.curselection()[0])
        name = w.get(index)
        for i, point in enumerate(self.points):
            if name == point.name:
                self.infobox.delete(0,END)
                self.infobox.insert(END,"NAME: " + point.name)
                self.infobox.insert(END,"PFAM: " + str(point.pfam.split(".")[0]))
                self.infobox.insert(END,"MOD: " + point.mod)
                self.infobox.insert(END,"NUM IN CLUSTER: " + str(int(self.sizes[i])))
                break

    def on_press(self, event):
        if self.toolbar._active is None:
            self.isPressed = True
            self.x0 = event.xdata
            self.y0 = event.ydata

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
                    if not self.points[i] in self.selectedpoints:
                        self.selectedpoints.append(self.points[i])
                        self.listbox.insert(END,str(self.points[i].name))

        self.rect.set_width(0)
        self.rect.set_height(0)
        self.rect.set_xy((self.x0, self.y0))
        try:
            self.scatter.figure.canvas.draw()
        except TypeError:
            pass


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
                        if not self.points[i] in self.selectedpoints:
                            self.selectedpoints.append(self.points[i])
                            self.listbox.insert(END,str(self.points[i].name))


    def on_key(self,event):
        if event.key == 'escape':
            self.selectedpoints = []
            self.selectedseqs = []
            self.listbox.delete(0, END)

        if event.key == 'r':
            print "Showing all selected points"

        if event.key =='n':
            self.scatter.clear()
            self.cmapindex += 1
            if self.cmapindex == len(self.cmaps):
                self.cmapindex = 0
            self.scatter.scatter(self.x_points, self.y_points,c=self.colors,cmap=self.cmaps[self.cmapindex],s=self.sizes)
            self.scatter.figure.canvas.draw()

        if event.key == 'c':
            self.scatter.clear()
            new_x = []
            new_y = []
            new_colors = []
            for (i, name) in enumerate(self.points):
                if self.points[i] in self.selectedpoints:
                    new_x.append(self.x_points[i])
                    new_y.append(self.y_points[i])
                    new_colors.append(self.colors[i])

            self.scatter.scatter(new_x, new_y,c=new_colors,s=self.sizes)
            self.scatter.add_patch(self.rect)
            self.scatter.figure.canvas.draw()

        if event.key == 'a':
            self.scatter.clear()
            self.scatter.scatter(self.x_points, self.y_points,c=self.colors,cmap=self.cmaps[0],s=self.sizes)
            self.scatter.add_patch(self.rect)
            self.scatter.figure.canvas.draw()

        if event.key == 's':
            with open(self.path, "w") as exportfile:
                for i, name in enumerate(self.selectedpoints):
                    exportfile.write(self.selectedpoints[i].line + '\n')
                    exportfile.write(self.selectedpoints[i].seq + '\n')
            print "Saved to file!"

    def connect(self):
        self.scatter.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.scatter.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.scatter.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.scatter.figure.canvas.mpl_connect('pick_event', self.on_pick)
        self.scatter.figure.canvas.mpl_connect('key_press_event', self.on_key)
        self.scatter.figure.canvas.mpl_connect('button_press_event',
                                               lambda event:self.scatter.figure.canvas._tkcanvas.focus_set())
        self.listbox.bind('<<ListboxSelect>>', self.onselect)



def tk_window_init(x_points,y_points):

    ## Create Tk main window and assign title and size
    root = Tk()
    root.wm_title("Clustering analysis results")
    root.geometry("1000x700+0+0")
    mainframe = Frame(master=root)
    mainframe.pack(fill=BOTH,expand=1)

    # Define the quit function used for the Quit Button.
    """# DO NOT USE X TO CLOSE, ONLY USE QUIT"""
    def _quit():
        root.quit()
        root.destroy()

    # Create matplotlib figure and main subplot (ax1)
    fig, ax1 = plt.subplots()
    # Plot points
    ax1.scatter(x_points, y_points,picker=2)
    # Remove axes and labels
    ax1.axes.get_xaxis().set_visible(False)
    ax1.axes.get_yaxis().set_visible(False)

    # Create the canvas and place matplotlib figure on canvas
    leftframe = Frame(master=mainframe)
    canvas = FigureCanvasTkAgg(fig, master=leftframe)
    canvas.show()

    #Create listbox and scrollbar
    rightframe = Frame(master=mainframe)
    listframe = Frame(master=rightframe)
    listbox = Listbox(master=listframe,selectmode=EXTENDED,height=30)
    listbox.pack(fill=BOTH,expand=1,anchor=CENTER,side=RIGHT)
    sb = Scrollbar(master=listframe,orient=VERTICAL)
    sb.configure(command=listbox.yview)
    listbox.configure(yscrollcommand=sb.set)
    listframe.pack(fill=BOTH,expand=1,side=TOP,padx=10)


    # Create listbox and infowindow for displaying selected proteins
    infobox = Listbox(master=rightframe,relief=SUNKEN,height=4)
    infobox.pack(fill=BOTH,expand=1,padx=10,pady=10)

    rightframe.pack(side=RIGHT,fill=BOTH,expand=1)
    leftframe.pack(side=LEFT,fill=BOTH,expand=1)

    # Import toolbar to navigate matplotlib
    toolbar = NavigationToolbar2TkAgg(canvas, leftframe)
    toolbar.update()

    # Create new frame for quit button
    bottomframe = Frame(master=rightframe,bg='blue')


    # Pack the components in a reasonable fashion
    toolbar.pack(side=TOP)
    bottomframe.pack(side=BOTTOM)
    canvas.get_tk_widget().pack(side=LEFT, fill=BOTH, expand=1)
    sb.pack(side=LEFT,fill=Y,pady=6,padx=3)

    # Create quit button
    quitbutton = Button(master=bottomframe, text='Quit', command=_quit)
    quitbutton.pack(side=LEFT,anchor=W)

    return root,ax1,toolbar,listbox,infobox
