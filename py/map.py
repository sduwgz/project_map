import sys

import init
import score

def DP(mole, ref):
    MIN_SCORE = -10000
    cols = len(ref)
    rows = len(mole)
    score_matrix = [[MIN_SCORE for col in range(cols + 1)] for row in range(rows + 1)]
    for i in xrange(len(score_matrix[0])):
        score_matrix[0][i] = 0
    map_matrix = [[-1 for col in range(cols + 1)] for row in range(rows + 1)]
    for i in xrange(rows):
        i += 1
        for j in xrange(cols):
            j += 1
            score_matrix[i][j] = score_matrix[i - 1][j - 1] + score.valid_score(mole[i - 1 : i], ref[j - 1 : j])
            map_matrix[i][j] = (i - 1) * cols + (j - 1)
            #DELETE: 3 is the max number of delete
            k = 0
            if j - 5 > 0:
                k = j - 5
            while k < j - 1:
                s = score_matrix[i - 1][k] + score.valid_score(mole[i - 1 : i], ref[k : j])
                if s > score_matrix[i][j]:
                    score_matrix[i][j] = s
                    map_matrix[i][j] = (i - 1) * cols + k
                k += 1
            #INSERT: 3 is the max number of insert
            k = 0
            if i - 5 > 0:
                k = i - 5
            while k < i - 1:
                s = score_matrix[k][j - 1] + score.valid_score(mole[k : i], ref[j - 1 : j])
                if s > score_matrix[i][j]:
                    score_matrix[i][j] = s
                    map_matrix[i][j] = k * cols + (j - 1)
                k += 1
    #the result is in the last row
    last_row = score_matrix[rows]
    best_score = last_row[1]
    end_site = 1
    #get the best score of last row
    for j in xrange(cols):
        j += 1
        if best_score < last_row[j]:
            best_score = last_row[j]
            end_site = j
    print '#####################'
    print end_site
    print '#####################'
 
    '''
    for i in xrange(rows):
        for j in xrange(cols):
            print '%f '%score_matrix[i][j],
        print
    
    #print the score_matrix
    for i in xrange(rows + 1):
        for j in xrange(cols + 1):
            print "%f "%score_matrix[i][j],
        print 
    '''

    #get the trace of DP
    #TODO
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print 'INPUT ERROR: Two files is needs. The first file is mole file(*.inx), and another is reference file(*.cmp)'
        exit()
    with file(sys.argv[1]) as f:
        moleset = init.read_moles(f)
    with file(sys.argv[2]) as f:
        ref = init.read_ref(f)[1 : ]
    #print 'ref: '
    #print ref
    for mole in moleset:
        print 'mole: %d'%mole
        #print moleset[mole]
        DP(moleset[mole], ref)
