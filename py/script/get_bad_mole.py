import sys
import string
from itertools import ifilter,imap

def get_bad_id(f):
    bad_id = []
    for line in ifilter(lambda x: len(x) > 0, imap(string.strip,f)):
        mole_id = abs(int(line.split()[0]))
        bad_id.append(mole_id)
    return bad_id

if __name__ == '__main__':
    bad_id = get_bad_id(sys.stdin)
    with file(sys.argv[1]) as f:
        line = f.readline().strip()
        while line:
            mole_id = int(line.split()[1])
            if mole_id in bad_id:
                print line
                line = f.readline().strip()
                print line
                line = f.readline().strip()
                print line
                line = f.readline().strip()
                print line
                line = f.readline().strip()
            else:
                line = f.readline().strip()
                line = f.readline().strip()
                line = f.readline().strip()
                line = f.readline().strip()
