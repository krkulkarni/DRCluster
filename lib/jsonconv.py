__author__ = 'kulkarnik'

def jsonmaker(colors,lines,matrix,fmat):

    colorset = set(colors)
    series = []
    for color in colorset:
        datalist = []
        parts = []
        for (i, line) in enumerate(lines):
            if (fmat == 'mod'):
                ## for modified format, namestrings are split at ";" character
                parts = line.split(";")
            elif (fmat == 'orig'):
                ## for original format, namestrings are split at "," character
                parts = line.split(",")

            if int(parts[1].strip()) == int(color):
                datalist.append(matrix.tolist()[i])

        jsonfile = {}
        jsonfile["name"] = color
        jsonfile["color"] = 'rgba(223, 83, 83, .5)'
        jsonfile["data"] = datalist
        series.append(jsonfile)


    return series