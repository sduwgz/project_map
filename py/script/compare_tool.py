import sys

pre_pos = 0
for line in sys.stdin:
    pos = float(line.split()[5])
    print pos - pre_pos
    pre_pos = pos
