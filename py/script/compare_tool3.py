import sys
import string

def get_score(f):
    d = {}
    for line in f.readlines():
        mole = line.split()[0]
        score = float(line.split()[1])
        d[mole] = score
    return d

if __name__ == '__main__':
    with file(sys.argv[1]) as f:
        d = get_score(f)
    for line in sys.stdin:
        line = line.strip()
        mole = line.split()[0]
        if d.has_key(mole) and float(line.split()[1]) > d[mole]:
            print line
