__author__ = 'kulkarnik'
import matplotlib.pyplot as plt
import time

def pyplotter2d(finalmat,mdstype,colors,namelist,t0):
    if (mdstype == 'm'):
        x_points = finalmat[:,0]
        y_points = finalmat[:,1]
    elif (mdstype == 'c'):
        x_points = []
        y_points = []
        for pair in finalmat:
            end = len(str(pair))
            x_points.append(str(pair)[2:end-2].strip().split()[0])
            y_points.append(str(pair)[2:end-2].strip().split()[1])


    fig, ax = plt.subplots()
    if (colors == []):
        colors = 'b'
    ax.scatter(x_points, y_points,c=colors)

    print len(namelist)
    print "Took", time.clock()-t0, "seconds"
    plt.show()
