def counter(a):
    count=0
    for i in a:
        if i.isalpha():
            count+=1
    return count
def constant_coord(aligned_seq,constant_loc):
    for i in range (constant_loc,len(aligned_seq)+1):
        if counter(aligned_seq[:i+1])==constant_loc+1 and aligned_seq[i]!='-':
            return i