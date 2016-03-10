#from sys inport argv
#script,imput_file = argv
import numpy as np

input_file = 'test_tmp'
current_file = open(input_file)
value = []
while 1:
    lines=current_file.readlines(100000)
    if not lines:
        break
    for line in lines:
        item = line.strip().split()
        if len(item)<=1:
            continue
        else:
            value.append(float(item[5])-float(item[4]))

print 'mean:',np.mean(value)
print 'variance:',(np.var(value))**0.5
print max(value),min(value)
print value
