import sys
from sklearn import cluster
sys.setrecursionlimit(100000)

def dbscan(matrix):
    db = cluster.DBSCAN(eps=1, min_samples=3).fit(matrix)
    labels = db.labels_
    return labels

class Object(object):
    def __init__(self,id,x,y):
        self.id = id
        self.x = x
        self.y = y

def _invert_dict_nonunique(d):
    newdict = {}
    for k, v in d.iteritems():
        newdict.setdefault(v, []).append(k)
    return newdict

def _find_points_near_point(coord,followed_coords,found,groups,all_coords,group,d,points):

    if coord.id in found:
        return found, groups, group, "f"

    followed_coords.add(coord.id)

    groups[points[coord.id]] = group    # group first point
    local_points = []

    for p in all_coords:
        ab_dist = ((coord.x - p.x)**2 +(coord.y - p.y)**2)**0.5
        if (ab_dist <= d):
            local_points.append(p)
    for p in local_points:
        if p.id == coord.id:
            continue

        if p.id not in followed_coords:
            found, groups, group, increase = _find_points_near_point(p,followed_coords,found,groups,all_coords,group,d,points)

    found.add(coord.id)
    return found, groups, group, "t"

def findgroups(matrix,d,points):
    indexnum = 0
    group = 1
    all_coords = []
    found = set()
    groups = {}
    groupids = []

    for line in matrix:
        o = Object(indexnum,float(line[0]),float(line[1]))
        indexnum += 1
        all_coords.append(o)

    for coord in all_coords:
        followed_coords = set()
        found, groups, group, increase = _find_points_near_point(coord,followed_coords,found,groups,all_coords,group,d,points)
        if (increase == "t"):
            group += 1
        groupids.append(groups[points[coord.id]])

    invgroups = _invert_dict_nonunique(groups)
    representatives, sizes = _findreps(invgroups)

    return groupids, representatives, sizes

def _findreps(invgroups):

    reps = []
    sizes = []

    print len(invgroups.keys()), "groups found"
    for key in invgroups.keys():
        found = False
        keypoints = invgroups[key]
        sizes.append(len(keypoints))
        for elem in keypoints:
            if elem.mod == 'pdb':
                reps.append(elem)
                found = True
                break
        if not found:
            for elem in keypoints:
                if elem.mod == 'mod':
                    reps.append(elem)
                    found = True
                    break
        if not found:
            reps.append(keypoints[0])

    return reps, sizes

def savereps(path, reps, distance):
    with open(path, 'w') as f:
        f.write("# " + str(len(reps)) + " groups found at distance " + str(distance) + "\n")
        for point in reps:
            f.write(point.line + "\n")
            f.write(point.seq + "\n")