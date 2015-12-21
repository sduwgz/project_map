import sys
import string
from itertools import ifilter, imap

def read_map(f):
    map_result = {}
    for line in ifilter(lambda x: len(x) > 0, imap(string.strip, f)):
        number = int(line.split()[0])
        start_site = int(line.split()[1])
        end_site = int(line.split()[2])
        map_result[number] = (start_site, end_site) 
    return map_result

if __name__ == '__main__':
    with file(sys.argv[1]) as f:
        huada_map = read_map(f)
    total = 0
    match = 0
    for line in ifilter(lambda x: len(x) > 0, imap(string.strip, sys.stdin)):
        if huada_map.has_key(int(line.split()[0])):
            total += 1
            start_site = int(line.split()[1])
            end_site = int(line.split()[2])
            if abs(start_site - huada_map[int(line.split()[0])][0]) < 3 or abs(end_site - huada_map[int(line.split()[0])][1]) < 3:
                match += 1
            else:
                print "%s "%line,
                print huada_map[int(line.split()[0])]
    print 'Total: %d'%total
    print 'Match: %d'%match
