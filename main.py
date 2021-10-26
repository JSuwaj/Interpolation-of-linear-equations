import numpy as np
import matplotlib.pyplot as plt
import time
def norm(res):
    ret=0
    for i in range(len(res)):
        ret+=res[i]*res[i]
    ret= np.sqrt(ret)
    return ret


def Jacobi(A,b):
    x=[]
    minres=10**-9
    normres=1
    k=0
    xOld=[]
    start=time.perf_counter()
    for i in range(len(A)):
        xOld.append(1)
    while normres>minres:
        x=[]
        for i in range(len(A)):
            sum=0
            for j in range(i):
                sum+=A[i][j]*xOld[j]
            for j in range(i+1, len(A)):
                sum += A[i][j]*xOld[j]
            x.append((b[i]-sum)/A[i][i])
        res=np.subtract(np.matmul(A,x),b)
        xOld=x
        normres=norm(res)
        if normres>10**20:
            print("algorytm jacobiego sie nie zbiega")
            return 0
        k+=1
    print("jacobi czas")
    print(time.perf_counter()-start)
    print("jacobi liczba iteracji")
    print(k)
    return time.perf_counter()-start


def GausSiedel(A,b):
    x=[]
    minres=10**-9
    normres=1
    k=0
    xOld=[]
    start=time.perf_counter()
    for i in range(len(A)):
        xOld.append(1)
    while normres>minres:
        x=[]
        for i in range(len(A)):
            sum=0
            for j in range(i):
                sum+=A[i][j]*x[j]
            for j in range(i+1, len(A)):
                sum += A[i][j]*xOld[j]
            x.append((b[i]-sum)/A[i][i])
        res=np.subtract(np.matmul(A,x),b)
        xOld=x
        normres=norm(res)
        if normres > 10 ** 20:
            print("algorytm gausa-siedla sie nie zbiega")
            return 0
        k+=1
    print("gaus siedel czas")
    print(time.perf_counter()-start)
    print("gaus siedel liczba iteracji")
    print(k)
    return time.perf_counter()-start


def LU(A,b):
    U = [row[:] for row in A]
    L=[[0 for i in range(len(A))] for j in range(len(A))]
    for i in range(len(A)):
        L[i][i]=1
    start = time.perf_counter()
    for k in range(len(A)):
        for j in range(k+1,len(A)):
            L[j][k]=U[j][k]/U[k][k]
            for l in range(k,len(A)):
                U[j][l]=U[j][l]-L[j][k]*U[k][l]
    y=[]
    for i in range(0,len(A)):
        sum=0
        for j in range(i):
            sum+=L[i][j]*y[j]
        y.append((b[i]-sum)/L[i][i])
    x = [0 for i in range(len(A))]

    x[len(A)-1]=y[len(A)-1]/U[len(A)-1][len(A)-1]
    for i in range(len(A)-1,-1,-1):
        sum = 0
        for j in range(i+1,len(A)):
            sum += U[i][j] * x[j]
        x[i]=(y[i]-sum)/U[i][i]
    res = np.subtract(np.matmul(A, x), b)
    normres = norm(res)
    print("LU norma z residum")
    print(normres)
    print("LU czas")
    print(time.perf_counter()-start)
    return time.perf_counter()-start

def ZadAB():
    N=927
    f=9
    a1=14
    a2=-1
    A=[[0 for i in range(N)] for j in range(N)]
    for i in range(N):
        for j in range(N):
            if i == j:
                A[i][j]=a1
            if i== j+1 or i==j-1 or i== j+2 or i==j-2:
                A[i][j]=a2
    b=[]
    for i in range(N):
        b.append(np.sin(i * (f + 1)))
    Jacobi(A,b)
    GausSiedel(A,b)


def ZadCD():
    N = 927
    f = 9
    a1 = 3
    a2 = -1
    A = [[0 for i in range(N)] for j in range(N)]
    for i in range(N):
        for j in range(N):
            if i == j:
                A[i][j] = a1
            if i == j + 1 or i == j - 1 or i == j + 2 or i == j - 2:
                A[i][j] = a2
    b = []
    for i in range(N):
        b.append(np.sin(i * (f + 1)))
    Jacobi(A, b)
    GausSiedel(A, b)
    LU(A, b)

def ZadE():
    N=[100,500,1000,2000,3000]
    f=9
    a1=14
    a2=-1
    LUTime=[]
    JacobiTime=[]
    GausSiedelTime=[]
    for k in range(5):
        A=[[0 for i in range(N[k])] for j in range(N[k])]

        for i in range(N[k]):
            for j in range(N[k]):
                if i == j:
                    A[i][j]=a1
                if i== j+1 or i==j-1 or i== j+2 or i==j-2:
                    A[i][j]=a2
        b=[]
        for i in range(N[k]):
            b.append(np.sin(i * (f + 1)))
        JacobiTime.append(Jacobi(A,b))
        GausSiedelTime.append(GausSiedel(A,b))
        LUTime.append(LU(A,b))
    plt.plot(N,JacobiTime)
    plt.plot(N,GausSiedelTime)
    plt.plot(N,LUTime)
    plt.legend(["Jacobi", "Gaus-Siedel","LU"])
    plt.xlabel("Size of a matrix")
    plt.ylabel("Time[s]")
    plt.show()

ZadAB()
ZadCD()
ZadE()
