def check(f_filename, g_filename):
    f = open(f_filename, "r")
    g = open(g_filename, "r")

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
            #if '24' in g_filename and not '8' in s and not '7' in s and len(s)<=6:
            #    print(s)
    
    for i in range(3,9):
        print(str(i)+':'+str(cnt[i]), end='   ')
    print('   sum='+str(sum(cnt)))

check("unsat_s.txt", "n_21_smallest_hulls.txt")
check("unsat_s.txt", "n_22_smallest_hulls.txt")
check("unsat_s.txt", "n_23_smallest_hulls.txt")
check("unsat_s.txt", "n_24_smallest_hulls.txt")
check("unsat_s.txt", "n_25_smallest_hulls.txt")
check("unsat_s.txt", "n_26_smallest_hulls.txt")
check("unsat_s.txt", "n_27_smallest_hulls.txt")
check("unsat_s.txt", "n_28_smallest_hulls.txt")
check("unsat_s.txt", "n_29_smallest_hulls.txt")
check("unsat_s.txt", "n_30_smallest_hulls.txt")
