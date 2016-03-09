from matplotlib.pyplot import step, show
import sys
import string
from itertools import ifilter, imap
import init

with file(sys.argv[1]) as f:
    ref = init.read_ref(f)
with file(sys.argv[2]) as f:
    moles= init.read_moles(f)
mole = moles[abs(int(sys.argv[3]))]
if int(sys.argv[3]) < 0:
    mole.reverse()
mole_length = len(mole)
start_site = 0
pre_length = 0
count = 0
print mole

x = range(mole_length)
for line in ifilter(lambda x: len(x) > 0, imap(string.strip, sys.stdin)):
    match_start = int(line.split()[2])
    print match_start
    #match_end = int(line.split()[3])
    ref_1 = ref[match_start : match_start + mole_length - count]
    ref_1 = mole[0 : count] + ref_1
  
    print ref_1
 
    for i in xrange(mole_length - 1):
        ref_1[i + 1] = ref_1[i + 1] + ref_1[i]
    count += 1
    step(x, ref_1, 'r', linewidth = 1)

for i in xrange(mole_length - 1):
    mole[i + 1] = mole[i + 1] + mole[i]
step(x, mole, 'b', linewidth = 3)
show()

