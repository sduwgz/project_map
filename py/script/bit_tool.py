import sys
import init

with file('../../data/ref.cmap') as f:
    ref = init.read_ref(f)
with file('../../data/input.bnx') as f:
    moles = init.read_moles(f)
state = 0
for line in sys.stdin.readlines():
    if state == 0:
        mole_id = int(line.split()[0])
        if mole_id < 0:
            state = 1
            continue
        start_site = int(line.split()[1])
        end_site = int(line.split()[2])
        print 
        print 'mole_id:%d'%mole_id
        mole = []
        for i in moles[mole_id]:
            mole.append(i / 100)
        print mole
        result = []
        print 'ref_nano:'
        for i in xrange(end_site - start_site + 1):
            result.append(ref[start_site + i] / 100)
        print result
        state = 1
        continue
    else:
        if mole_id < 0:
            state = 0
            continue
        start_site = int(line[1:-2].split(',')[0])
        end_site = int(line[1:-2].split(',')[1])
        result = []
        print 'ref_huada:'
        for i in xrange(end_site - start_site + 1):
            result.append(ref[start_site + i] / 100)
        print result
        state = 0
        continue
