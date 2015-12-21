import sys
import init

with file('../../data/ref.cmap') as f:
    ref = init.read_ref(f)
with file('../../data/input.bnx') as f:
    moles = init.read_moles(f)
mole_id = int(sys.argv[1])
print 'mole:'
print moles[mole_id]
start_site = int(sys.argv[2])
end_site = int(sys.argv[3])
result = []
print 'ref:'
for i in xrange(end_site - start_site + 1):
    result.append(ref[start_site + i])
print result
