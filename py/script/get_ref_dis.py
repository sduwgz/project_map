import init
import sys

if __name__ == '__main__':
    start = int(sys.argv[1])
    end = int(sys.argv[2])
    with file('/home/wei/wei/study/nano/nanoMap/data/ref_filter.cmap') as f:
        dis = init.read_ref(f)
    i = start
    while i < end:
        print '%d'%dis[i],
        i += 1
    
