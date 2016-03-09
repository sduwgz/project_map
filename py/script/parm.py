import sys

number = []
s = 0
for line in sys.stdin:
    a = int(line.split()[4])
    b = int(line.split()[5])
    if abs(a - b) < 1000:
        s += a
        number.append(b - a)
print (s + 0.0) / len(number)
s = 0
for i in number:
    s += i
mean = (s + 0.0) / len(number)
s = 0
for i in number:
    s += (i - mean)**2
var = s / len(number)
print mean
print var**(0.5)
