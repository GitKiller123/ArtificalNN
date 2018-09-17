import numpy as np
import time

start_time = time.time()

errlist = []
trialstot = 100000

n = 100
p = np.array([12,20,40,60,80,100])
for choice in p:
    errcount = 0
    currp = np.zeros([choice,n])
    for trials in range(trialstot):
        w = np.zeros([n, n])
        for x in range(choice):
            currp[x, :] = np.random.choice([-1, 1], n)
            w = w + np.outer(currp[x, :], currp[x, :])/n

        np.fill_diagonal(w, 0)
        pnum = np.random.randint(0,choice)
        nnum = np.random.randint(0,100)

        nout = w[nnum, :].dot(currp[pnum, :])

        if nout < 0:
            nout = -1
        else:
            nout = 1
        if nout == currp[pnum, nnum]:
            continue
        else:
            errcount = errcount + 1
    errperc = errcount / trialstot
    errlist.append(errperc)
    print(errcount)
print(errlist)

print("--- %s seconds ---" % (time.time() - start_time))