import random
import sys
def get_random_seq(l):
    seq = ''
    i = 0
    d = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}
    while i < l:
        a = random.randint(0, 3)
        seq += d[a]
        i += 1
    return seq

def get_mole_site(pattern, seq):
    k = 0
    while True:
        k = seq.find(pattern, k + 1)
        if k == -1:
            break
        print '1	4639675.0	683	1	1	%d	1.0	1	1'%k
if __name__ == '__main__':
    l = int(sys.argv[1])
    seq = get_random_seq(l)
    get_mole_site('GCCGTT', seq)
