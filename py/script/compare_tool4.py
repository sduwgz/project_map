import sys
import string

if __name__ == '__main__':
    good_map = []
    for line in sys.stdin:
        good_map.append(int(line.split()[0]))
    f = file(sys.argv[1])
    line = f.readline()
    while line:
        if line.startswith('m'):
            if int(line[4:]) in good_map:
                line = f.readline()
                while line:
                    if line.startswith('m'):
                        break
                    else:
                        print line,
                        line = f.readline()
            else:
                line = f.readline()
        else:
            line = f.readline()
