import sys
import string

def filter_ref(f):
    lines = f.readlines()
    i = 0
    count = 1
    while i < len(lines) - 3:
        if lines[i].startswith('#'):
            #print lines[i].strip()
            i += 1
            continue
        pos1 = float(lines[i].split()[5])
        pos2 = float(lines[i + 1].split()[5])
        pos3 = float(lines[i + 2].split()[5])
        pos4 = float(lines[i + 3].split()[5])
        if pos2 - pos1 > 600:
            print '1   4639675.0   683 %d 1   %f   1.0 1   1'%(count, float(pos1))
            i += 1
            count += 1
        if pos2 - pos1 < 600 and pos3 - pos2 > 600:
            print '1   4639675.0   683 %d 1   %f   1.0 1   1'%(count, float(int((pos1 + pos2) / 2 + 0.5) ) )
            i += 2
            count += 1
        if pos2 - pos1 < 600 and pos3 - pos2 < 600 and pos4 - pos3 > 600:
            print '1   4639675.0   683 %d 1   %f   1.0 1   1'%(count, float(int((pos1 + pos2 + pos3) / 3 + 0.5) ) )
            i += 3
            count += 1
        if pos2 - pos1 < 600 and pos3 - pos2 < 600 and pos4 - pos3 < 600:
            print '1   4639675.0   683 %d 1   %f   1.0 1   1'%(count, float(int((pos1 + pos2 + pos3 + pos4) / 4 + 0.5) ) )
            i += 4
            count += 1

if __name__ == '__main__':
    with file(sys.argv[1]) as f:
        filter_ref(f)
#    pos = []
#    for line in imap(string.strip,isys.stdin).readlines():
#        if line.startswith('#'):
#            continue
#        pos.append(float(line.split()[5]))
#    for 
