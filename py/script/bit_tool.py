import sys
import init

with file(sys.argv[1]) as f:
    ref = init.read_ref(f)
with file(sys.argv[2]) as f:
    moles = init.read_moles(f)
state = 0
for line in sys.stdin.readlines():
    if state == 0:
        mole_id = int(line.split()[0])
        start_site = int(line.split()[2])
        end_site = int(line.split()[3])
        print 
        print 'mole_id:%d'%mole_id
        mole = []
        for i in moles[abs(mole_id)]:
            mole.append(i)
        print mole
        result = []
        print 'ref_nano:'
        for i in xrange(end_site - start_site + 1):
            result.append(ref[start_site + i])
        if mole_id < 0:
            result.reverse()
        print result
        state = 1
    else:
        start_site = int(line[1:-2].split(',')[0])
        end_site = int(line[1:-2].split(',')[1])
        result = []
        print 'ref_huada:'
        for i in xrange(end_site - start_site + 1):
            result.append(ref[start_site + i])
        print result
        state = 0
