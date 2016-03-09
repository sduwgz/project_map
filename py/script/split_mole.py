import sys

import init

if __name__ == '__main__':
    with file(sys.argv[1]) as f:
        moles = init.read_moles(f)
    mole = moles[int(sys.argv[2])]
    start = int(sys.argv[3])
    for i in xrange(len(mole) - 4):
        print '0 %d'%int(sys.argv[2])
        site = start
        print '0',
        for j in xrange(4):
            print '\t%d'%site,
            site += mole[i + j]
        print 
        print 'QX11'
        print 'QX12'
        start += mole[i]
