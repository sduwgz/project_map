import sys
import string
import init

with file(sys.argv[1]) as f:
    ref = init.read_ref(f)
with file(sys.argv[2]) as f:
    mole_set = init.read_moles(f)

line = sys.stdin.readline()
line = line.split()

mole_id = int(line[1])
mole_direct = line[2]
mole = mole_set[mole_id]
print mole
pos = 18
d = int(sys.argv[3])
while pos < len(line):
    start = int(line[pos])
    left = int(line[pos + 1])
    right = int(line[pos + 2])
    if left == -1 or right == -1:
        pos += 3
        continue
    else:
        pos += 3
        break

while pos < len(line):
    start_next = int(line[pos])
    left_next = int(line[pos + 1])
    right_next = int(line[pos + 2])
    if left_next == -1 or right_next == -1:
        pos += 3
        continue
    mole_length = 0
    ref_length = 0
    i = start - 1
    while i < start_next - 1:
        mole_length += mole[i]
        i += 1
    j = left + d
    while j < right_next + d:
        ref_length += ref[j]
        j += 1
    print "%d %d %d %d %d %d"%(start, start_next - 1, left, left_next - 1, mole_length, ref_length)
    start = start_next
    left = right_next
    pos += 3
