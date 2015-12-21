import sys

for line in sys.stdin:
    l = line.split()
    pre_pos = 0
    for i in l:
        print float(i) - pre_pos
        pre_pos = float(i)
