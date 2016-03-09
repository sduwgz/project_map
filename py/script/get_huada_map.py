import sys
import string
from itertools import ifilter, imap

def pos_to_site(f):
    d = {}
    for line in f.readlines():
        if line.startswith('#'):
            continue
        site = int(float(line.split()[3]))
        pos = int(float(line.split()[5]))
        d[pos] = site
    return d

if __name__ == '__main__':
    with file(sys.argv[1]) as f:
        d = pos_to_site(f)
    for line in ifilter(lambda x: len(x) > 0, imap(string.strip, sys.stdin)):
        mole_id = line.split()[1]
        de = line.split()[6]
        start_site = int(float(line.split()[9]))
        end_site = int(float(line.split()[10]))
        if d.has_key(start_site) and d.has_key(end_site):
            print '%s %d %d'%(de + mole_id, d[start_site], d[end_site])
