import math
import qsharp
import json
import timeit

from Shor import FindOrder

N = 33
X = 10

N_BIT = math.ceil(math.log(N+1, 2)) 
m = []

tic = timeit.default_timer()

# for i in range(1000):
#     m.append(''.join(str(x) for x in reversed(
#         FindOrder.simulate(N = N, n = N_BIT, Nr = X)
#     )))

print(FindOrder.simulate(N = N, n = N_BIT, Nr = X, shots = 1))

toc = timeit.default_timer()

print(toc - tic)

# res = {}
# for i in range(2**(2*N_BIT)):
#     res[i] = 0

# for m0 in m:
#     m0  = int(m0,2)
#     res[m0]+= 1

# print(res)
