import os


def tparse(file_name, to_report, dlim='\t', skip='#', sf=False, fter=None):
    """ 
    parse table file
    to_report is a list [2,4,6] (0-based)
    sf : skip first line
    """
    if type(to_report) is int:
        REPORTONE=True
    else:
        REPORTONE=False

    with open(file_name) as f:
        if sf:
            f.readline()

        for line in f:
            if line.startswith(skip):
                continue

            fields = line.strip().split(dlim)

            if fter and fter(fields):
                continue
                
            if REPORTONE:
                yield fields[to_report]
            else:
                yield [fields[ind] for ind in to_report]


def ldir(dir_name, ext):
    """ 
    list dir
    ext is the extension of the file
    """

    for f in os.listdir(dir_name):
    	if f.endswith(ext):
	    yield f

def wdir(dir_name, ext):
    """ 
    dirnames are the names of the subfolders
    filenames are the names of the files in the current folder
    """
    for dirpath, dirnames, filenames in os.walk(dir_name):
        for valid_file in filenames:
            if valid_file.endswith(ext):
                yield os.path.join(dirpath, valid_file)


def registerlist(d, key, value):
    if key in d:
        d[key].append(value)
    else:
        d[key] = [value]


def registercnt(d, key):
    if key in d:
        d[key] += 1
    else:
        d[key] = 1

    return

def registerset(d, key, value):
    if key in d:
        d[key].add(value)
    else:
        d[key] = set([value])

def increment_cnt(d, key, cnt):
    if key in d:
        d[key] += cnt
    else:
        d[key] = cnt

def registerkey_dict(d, key):
    if key not in d:
        d[key] = {}

    return


def pathloc(path, loc):

    """
    return the location of a path
    """
    return path.split('/')[loc]
