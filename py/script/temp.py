import sys

import init

if __name__ == '__main__':
    with file(sys.argv[1]) as f:
        moles = init.read_moles(f)
    s = 0
    count = 0
    for i in moles:
        for j in moles[i]:
            s += j
            count += 1
    print s / count
