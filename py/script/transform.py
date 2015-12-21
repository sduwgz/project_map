import sys

with file(sys.argv[1]) as f:
    line = f.readline()
    while line:
        number = int(line.split()[0])
        score = float((line.split()[1]))
        start_site = int(line.split()[2])
        end_site = int(line.split()[3])
        line = f.readline()
        if float(line.split()[1]) > score:
            number = int(line.split()[0])
            start_site = int(line.split()[2])
            end_site = int(line.split()[3])
        if number > 0:
            print '+%d '%number,
            print ' %d %d'%(start_site, end_site)
        else:
            print '%d '%number,
            print ' %d %d'%(start_site, end_site)
        line = f.readline()
