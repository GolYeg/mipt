import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg
import time


def MakeCoefMatr(h, c):
    depth = int (1 / h)
    netSize = depth + 1
    eqNum = (depth - 1) * (depth - 1)
    coefMatr = np.zeros((eqNum, eqNum), dtype='float64')

    position = 0
    for j in range(1, netSize - 1):
        for i in range(1, netSize - 1):
            coefMatr[position, position] = 4 + c * (h ** 2)
            if i != depth - 1:
                coefMatr[position, position + 1] = -1
            if i != 1:
                coefMatr[position, position - 1] = -1
            if j != depth - 1:
                coefMatr[position, position + (depth - 1)] = -1
            if j != 1:
                coefMatr[position, position - (depth - 1)] = -1
            position += 1
    #print (zeroMatr)

    rPart = np.empty(eqNum)
    rPart.fill(h ** 2)

    return coefMatr, rPart, netSize

def MakeLuView(coefMatr, rPart, netSize):
    start = time.time()

    lu, piv = linalg.lu_factor(coefMatr)
    x = linalg.lu_solve((lu, piv), rPart)

    finish = time.time()
    #print('time - ', finish - start)
 
    matrOfVal = np.zeros((netSize, netSize), dtype='float64')
    position = 0
    for j in range(1, netSize - 1):
        for i in range(1, netSize - 1):
            matrOfVal[i, j] = x[position]
            position += 1

    return matrOfVal

def TimeCheck():
    flag = 0
    hParams = []
    times = []
    backH = 4
    while flag == 0:
        start = time.time()
        coefMatr, rPart, netSize = MakeCoefMatr(1 / backH, 0.1)
        MakeLuView(coefMatr, rPart, netSize)
        finish = time.time()

        hParams.append(backH)
        times.append(finish - start)
        if backH < 200:
            backH += 10
        else:
            break
        #print(backH + 1)
        #print(finish - start)

    plt.plot(hParams, times)
    plt.xlabel('Size of net')
    plt.ylabel('Time of calculation')
    plt.grid()
    plt.show()

TimeCheck() 

start = time.time()
coefMatr, rPart, netSize = MakeCoefMatr(1/200, 0.1)

# plt.spy(coefMatr, markersize=5, marker='.', color='blue')
# plt.grid()
# plt.title('Coefficient Matrix')
# plt.show()

# #print('CoefMatr - ', coefMatr, '\n')
# #print('RPart - ', rPart, '\n')

matrOfVal = MakeLuView(coefMatr, rPart, netSize)
finish = time.time()
print (finish - start)
#print (matrOfVal)
plt.imshow(matrOfVal, cmap='rainbow')
plt.savefig ("rainbow diogram")     


