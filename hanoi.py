"""
Busybar widget
Rick Lawson
r_b_lawson at yahoo dot com
Heavily borrowed from ProgressBar.py which I got off the net but can't remember where
Feel free to add credits.
Comments by Stewart Midwinter: 
 I added a Quit button so you can stop the app.
 I also set up a timer so that the BusyBar stops after a certain period.
 Next, I added a button to bring up a BusyBar in a top-level window, similar to what
 you might a process to do while it was, in fact, busy.
 The top-level window is non-modal; it's left as an exercise for you, the reader, to change that if needed.


config options
--------------
BusyBar is derived from frame so all frame options are fine
Here are the options specific to this widget
fill       - color of the progress box (the box that bounces back and forth)
boxWidth   - width of progress box as a fraction of total widget width
interval   - interval in ms at which the progress box is moved
             ie, the shorter this is the faster the box will move
increment  - fraction of widget width that the box moves during an update
             ie, 0.05 means that box will move 5% of the total width at an update
text       - text of message that is displayed in the middle of the widget
foreground - color of text message
font       - font of text message
"""

from Tkinter import *
import time

class BusyBar(Frame):
    def __init__(self, master=None, **options):  
        # make sure we have sane defaults
        self.master=master
        self.options=options
        self.width=options.setdefault('width', 100)
        self.height=options.setdefault('height', 10)
        self.background=options.setdefault('background', 'gray')
        self.relief=options.setdefault('relief', 'sunken')
        self.bd=options.setdefault('bd', 2)
        
        #extract options not applicable to frames
        self._extractOptions(options)
        
        # init the base class
        Frame.__init__(self, master, options)
        
        self.incr=self.width*self.increment
        self.busy=0
        self.dir='right'
        
        # create the canvas which is the container for the bar
        self.canvas=Canvas(self, height=self.height, width=self.width, bd=0,
                           highlightthickness=0, background=self.background)
        # catch canvas resizes
        self.canvas.bind('<Configure>', self.onSize)
        
        # this is the bar that moves back and forth on the canvas
        self.scale=self.canvas.create_rectangle(0, 0, self.width*self.barWidth, self.height, fill=self.fill)
                                                
        # label that is in the center of the widget
        self.label=self.canvas.create_text(self.canvas.winfo_reqwidth() / 2,
                                           self.height / 2, text=self.text,
                                           anchor="c", fill=self.foreground,
                                           font=self.font)
        self.update()
        self.canvas.pack(side=TOP, fill=X, expand=NO)
        
    def _extractOptions(self, options):
        # these are the options not applicable to a frame
        self.foreground=pop(options, 'foreground', 'yellow')
        self.fill=pop(options, 'fill', 'blue')
        self.interval=pop(options, 'interval', 30)
        self.font=pop(options, 'font','helvetica 10')
        self.text=pop(options, 'text', '')
        self.barWidth=pop(options, 'barWidth', 0.2)
        self.increment=pop(options, 'increment', 0.05)

    # todo - need to implement config, cget, __setitem__, __getitem__ so it's more like a reg widget
    # as it is now, you get a chance to set stuff at the constructor but not after
        
    def onSize(self, e=None):
        self.width = e.width
        self.height = e.height
        # make sure the label is centered
        self.canvas.delete(self.label)
        self.label=self.canvas.create_text(self.width / 2, self.height / 2, text=self.text,
                                           anchor="c", fill=self.foreground, font=self.font)

    def on(self):
        self.busy = 1
        self.canvas.after(self.interval, self.update)
        
    def of(self):
        self.busy = 0

    def update(self):
        # do the move
        x1,y1,x2,y2 = self.canvas.coords(self.scale)
        if x2>=self.width:
            self.dir='left'
        if x1<=0:
            self.dir='right'
        if self.dir=='right':
            self.canvas.move(self.scale, self.incr, 0)
        else:
            self.canvas.move(self.scale, -1*self.incr, 0)

        if self.busy:
            self.canvas.after(self.interval, self.update)
        self.canvas.update_idletasks()
        
def pop(dict, key, default):
    value = dict.get(key, default)
    if dict.has_key(key):
        del dict[key]
    return value
        
        
if __name__=='__main__':
    root = Tk()
    
    def popup():
        win=Toplevel()
        win.title("I'm busy too!")
        bb1=BusyBar(win, text='Wait for me!')
        bb1.pack()
        for i in range(0,30):
                time.sleep(0.1)
                bb1.update()
                root.update()
        bb1.of()
        time.sleep(1)
        win.destroy()

    t = Text(root)
    t.pack(side=TOP)
    bb = BusyBar(root, text='Please Wait')
    bb.pack(side=LEFT, expand=NO)
    but = Button(root, text= 'Pop-up BusyBar', command=popup)
    but.pack(side=LEFT, expand=NO)
    q = Button(root, text= 'Quit', command=root.destroy)
    q.pack(side=LEFT, expand=NO)
    l = Label(root, text="I'm a status bar !")
    l.pack(side=RIGHT)
    bb.on()
    root.update_idletasks()
    for i in range(0,30):
        time.sleep(0.1)
        root.update()
    bb.of()
    root.mainloop()