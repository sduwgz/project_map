import sys
import string
from itertools import ifilter, imap

def read_ref(f):
    ref = {}
    for line in ifilter(lambda x: len(x) > 0, imap(string.strip, f)):
        if line.startswith('#'):
            continue
        else:
            number = int(line.split()[3])
            pos = int(float(line.split()[5]))
            ref[pos] = number
    return ref

def read_map(f):
    map_label = []
    mapp = []
    for line in ifilter(lambda x: len(x) > 0, imap(string.strip, f)):
        mole_id = line.split()[1]
        if int(line.split()[11]) == int(line.split()[12]) + int(line.split()[13]): 
            direction = line.split()[6]
            begin_pos = int(line.split()[9])
            end_pos = int(line.split()[10])
            key = direction + mole_id
            mapp.append((begin_pos, end_pos))
            map_label.append(key)
    return mapp, map_label

if __name__ == '__main__':
    with file(sys.argv[1]) as f:
        ref = read_ref(f)
    with file(sys.argv[2]) as f:
        mapp, map_label = read_map(f)
    for i in xrange(len(mapp)):
        start_pos = mapp[i][0]
        end_pos = mapp[i][1]
        if ref.has_key(start_pos):
            start_site = ref[start_pos]
        else:
            continue
        if ref.has_key(end_pos):
            end_site = ref[end_pos]
        else:
            continue
        print map_label[i],
        print " %d  %d"%(start_site, end_site)
