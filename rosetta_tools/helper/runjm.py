import sys

class InfoObject():
    qlen = 1
    domainnum = 0
    alnstart = 0
    alnend = 0
    score = 0
    alnseq = ''
    description = None

def parsedom(domhandle,dict):

    try:
        while True:
            domline = domhandle.next().split()
            if (domline[0].startswith('#')):
                continue
            targetname = domline[0]
            newobj = InfoObject()
            if not targetname in dict:
                dict[targetname] = []
            else:
                pass

            newobj.domainnum = int(domline[9])
            newobj.qlen = int(domline[5])
            newobj.alnstart = int(domline[17])
            newobj.alnend = int(domline[18])
            newobj.description = domline[22]
            newobj.score = float(domline[13])

            dict[targetname].append(newobj)

    except StopIteration:
        print("done parsing domfile")

    return dict

def parsealn(alnhandle,dict):

    try:
        while True:
            alnline = alnhandle.next()
            if (alnline.startswith('#') or alnline == '\n' or alnline == '//\n'):
                continue
            splitline = alnline.split()
            name = splitline[0].split('/')[0]
            try:
                range = splitline[0].split('/')[1]
                rangestart = int(range.split('-')[0])
                rangeend = int(range.split('-')[1])

            except IndexError:
                newobj = InfoObject()
                dict[name] = []
                newobj.alnstart = 1
                newobj.alnend = 1
                rangestart = 1
                rangeend = 1
                dict[name].append(newobj)

                # add an entry to dict for query name

            templist = dict.get(name)
            for obj in templist:
                if (obj.alnstart == rangestart and obj.alnend == rangeend):
                    obj.alnseq = splitline[1]
                    break

    except StopIteration:
        print("done parsing alnfile")

def cleanstring(fullstring):
    return fullstring.replace('-','')

def objinfo(dict,pdblog):

    seqlist = []
    cleanseq = []
    keys = dict.keys()
    for key in keys:
        templist = dict.get(key)
        for obj in templist:
            clean = cleanstring(obj.alnseq)

            if (((obj.score/obj.qlen) > 1.3) and clean not in cleanseq):
                temp = str(obj.alnstart) + ' ' + str(obj.alnseq)
                accanddomain = key.split(":")
                print accanddomain[0]+"_"+accanddomain[1].split(",")[0],
                pdblog.write(accanddomain[0] + ' ' + accanddomain[1].split(",")[0]+"\n")
                cleanseq.append(clean)
                seqlist.append(temp)

    return seqlist

def runjm(fastafile,domfile,alignfile):
    seqhandle = open(fastafile, 'rt')
    domhandle = open(domfile,'rt')
    alnhandle = open(alignfile, 'rt')
    infodict = {}

    orig_std_out = sys.stdout

    parsedom(domhandle,infodict)
    parsealn(alnhandle,infodict)

    with open('outlog', 'w') as log:
        with open('pdblog', 'w') as pdblog:
            sys.stdout = log
            lines = seqhandle.readlines()
            seqlist = []
            print '##',
            for line in lines:
                if (line[0] == '>'):
                    #print line[1:].strip(),
                    print sys.argv[1].split(".")[0],
                    seqlist.append(str(1) + " " + infodict[line[1:].strip()][0].alnseq)
            seqlist.extend(objinfo(infodict,pdblog))

            print ''
            for item in seqlist:
                print item



    sys.stdout = orig_std_out
    seqhandle.close()
    domhandle.close()
    alnhandle.close()

if (__name__ == "__main__"):
    runjm(sys.argv[1],sys.argv[2],sys.argv[3])
    #      fasta        hits        align
