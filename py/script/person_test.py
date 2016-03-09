import scipy.stats as st

v = [564, 2849, 1618, 856, 259, 129, 83, 27, 76]

p = []
t = 6458
for i in range(8):
    p.append(st.poisson.pmf(i, 2) * 6458)
    t -= st.poisson.pmf(i, 2) * 6458
p.append(t)
sum = 0
print p
for i in range(9):
    sum += v[i]**2 / p[i]
print sum - 6458



