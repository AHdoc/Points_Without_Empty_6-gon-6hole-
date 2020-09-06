def sz(s):
    rt = 0
    x = int(s)
    while x != 0:
        rt += x % 10
        x = x // 10
    return rt


f = open("unsat_s.txt", "r")
unsat_s = [[] for i in range(31)]
for line in f.readlines():
    for s in line.split():
        unsat_s[sz(s)].append(s)
f.close()


g = open("unsat_s.txt", "w")
for i in range(1,31):
    unsat_s[i] = list(set(unsat_s[i]))
    unsat_s[i].sort(reverse=True)
    for x in unsat_s[i]:
        g.write(x+' ')
    g.write("\n")
g.close()
