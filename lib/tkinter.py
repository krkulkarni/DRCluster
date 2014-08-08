__author__ = 'kulkarnik'

from Tkinter import Listbox,Tk,END
import tkMessageBox

top = Tk()

Lb1 = Listbox(top)
Lb1.insert(END, "Python")

Lb1.pack()
top.mainloop()