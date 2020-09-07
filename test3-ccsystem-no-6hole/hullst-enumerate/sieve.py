f = open("unsat_s.txt", "r")
g = open("n_30_smallest_hulls.txt", "r")

unsat_s = []
for line in f.readlines():
    for s in line.split():
        unsat_s.append(s)

cnt = [0 for i in range(10)]
for line in g.readlines():
    tmp = line.split()
    s = tmp[2]

    chk = True
    for s2 in unsat_s:
        if s2 in s:
            chk = False
    
    if chk:
        cnt[ord(s[0])-ord('0')] += 1
        #print(s)

for i in range(3,9):
    print(str(i)+':'+str(cnt[i]), end='   ')
print()
print(sum(cnt))
