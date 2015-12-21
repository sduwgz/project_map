import math
import scipy.stats as stats

def get_score(mole_blocks, ref_blocks):
    assert(len(mole_blocks) == len(ref_blocks))
    score = 0.0
    for i in xrange(len(mole_blocks)):
        mb = mole_blocks[i]
        rb = ref_blocks[i]
        score += valid_score(mb, rb)
    return score

def valid_score(mole_block, ref_block):
    mb = mole_block
    rb = ref_block
    mb_length = 0
    rb_length = 0
    miss_site = 0
    for i in mb:
        mb_length += i
    for i in rb:
        rb_length += i
        if i < 1000:
            miss_site += 1
    delta = mb_length - rb_length
    if len(mb) > 1:
        #log(pI(number)) is the punish of score
        score = guss(delta) + pI(len(mb) ) - background(delta)
    elif len(rb) - miss_site > 1:
        delete = int(((len(rb) - miss_site - 1 + 0.0) / mb_length * 10000 ) + 0.5)
        if delete < 1:
            delete = 1
        if delete > 20:
            delete = 20
        #log(pD(delete)) is the punish of score
        score = guss(delta) + pD(delete) - background(delta)
    else:
        # log(1) means no punish 
        score = guss(delta) + math.log(1) - background(delta)
    return score
    
def guss(delta):
    mu = 46.9
    sigma = 570.0
    return 0.0 - (0.5 * math.log(2) +0.5 * math.log(math.pi) + math.log(sigma) ) - (delta - mu) ** 2 / (2 * sigma ** 2)
    #return 0.0 - (math.log((2 * math.pi) ** 0.5 * sigma) + (delta - mu) ** 2 / (2 * sigma ** 2) )

def pI(k):
    lambd = 1.064
    #return stats.expon.pdf(k, lambd)
    #return lambd * math.exp(-lambd * k)
    return math.log(lambd) - lambd * k

def pD(k):
    lambd = 1.82
    #return stats.poisson.pmf(k, lambd)
    return k * math.log(lambd) - lambd - math.log(fact(k) )

def fact(k):
    f = 1
    for i in xrange(k):
        f *= (i + 1)
    return f


def background(delta):
    mu = 1870.0
    sigma = 10840
    return 0.0 - (0.5 * math.log(2) +0.5 * math.log(math.pi) + math.log(sigma) ) - (delta - mu) ** 2 / (2 * sigma ** 2)
    #return 0.0 - (math.log((2 * math.pi) ** 0.5 * sigma) + (delta - mu) ** 2 / (2 * sigma ** 2) )


if __name__ == '__main__':
    mole = [[1825,1000],[1743],[8081],[3836],[3396],[5158],[2698],[3864],[5521],[5968],[10000,6729],[13970],[21305, 10000],[12926]]
    ref = [[2837],[2062],[7963],[4015],[3526],[5188],[2833],[3831],[5280],[6240],[16985],[10000,4055],[31788],[12825]]
    r1 = [[5900,800],[14110],[7600,16400],[1230,470,6200],[2330,620,12400],[3000,4700],[7300],[14100]]
    r2 = [[6700],[14400],[24800],[750,10900],[15000],[8000],[18400],[3340]]
    m1 = [[6934],[14541],[25342],[11226],[14923],[8040],[6489],[10582,3424]]
    m2 = [[6934],[14541],[25342],[11226],[14923],[8040],[17282],[3424]]
    
    #print get_score(m1, r1)
    #print get_score(m2, r2)
    m = [8000]
    r = [7000, 500]
    print valid_score(m, r)
    print pD(5)
