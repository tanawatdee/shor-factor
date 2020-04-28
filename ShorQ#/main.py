import math
import qsharp
import json
import random
import timeit

from Shor import FindOrder

def gcd(a, b):
    while b:
        a, b = b, a%b
    return a

def cont_frac(k, r):
    seq = []
    while r:
        seq.append(k//r)
        k,r = r, k%r
    return seq

def simp_frac(seq):
    k,r = 0,1
    for i in reversed(seq):
        k,r = r, k+r*i
    return r,k

def mod_inverse(m, a):   
    x,y,z = 1,0,m
    while a > 1:  
        q,m,a = a//m, a%m, m
        x,y = y, x-q*y
    return x+z if x<0 else x
            
def factor(N):                          # number to factor
    N_BIT = math.ceil(math.log(N+1, 2)) # N is represented in N_BIT bits
    N_QFT = 2**(2*N_BIT)                # number of basis states in QFT
    p,q = None,None                     # 1st factor, 2nd factor
    r = {                               # return value
        'factor': None,
        'log'   : []
    }

    while p is None or q is None:

        if N%2 == 0: # quantum part is not designed for even number
            p,q = 2, N//2
            r['log'].append({ 'type': 'DIV2' })
            r['factor'] = [p,q]
            print(N, 'is divisible by 2.')
            break

        X = random.randint(2, N-1) # chooe a random number 1<X<N

        print('---------------------------')
        print('choose X =', X)

        gcdXN = gcd(X,N)

        if gcdXN != 1:
            p,q = gcdXN, N//gcdXN
            r['log'].append({
                'type': 'GCD',
                'rand': X
            })
            r['factor'] = [p if p < q else q, q if q > p else p]
            print(X, 'and', N, 'are not coprime.')
            print('apply gcd(', X, ',', N, ')')
            break
            
        log = {
            'type': 'QFT',
            'rand': X,
            'measure': []
        }

        ############## This is quantum part #############

        m = []

        tic = timeit.default_timer()

        for i in range(2*N_BIT):
            m.append(''.join(str(x) for x in reversed(
                FindOrder.simulate(N = N, n = N_BIT, Nr = X)
            )))

        toc = timeit.default_timer()
        r['time'] = toc - tic

        #################### end of quantum part ##################

        for m0 in m:
            m0  = int(m0,2)
            seq = cont_frac(m0, N_QFT)
            log_m0 = {
                'value': m0,
                'basis': N_QFT,
                'frac' : []
            }
            print('')
            print('measured value:', m0, '(of', N_QFT, 'basis states)')
            print('continued fraction:', seq)

            for j in range(1, len(seq)):
                k0,r0 = simp_frac(seq[:j+1])

                if r0 > N:
                    break

                if abs(k0/r0 - m0/N_QFT) > 1/N_QFT:
                    continue

                log_m0['frac'].append({ 'k': k0, 'r': r0 })
                print('guessed k/r:', k0, '/', r0)

                if r0%2 or X**(r0//2)%N == -1:
                    continue

                Xpower = X**(r0//2)
                Xplus  = Xpower + 1
                Xminus = Xpower - 1
                p_plus  = gcd(Xplus, N)
                p_minus = gcd(Xminus, N)
                q_plus  = N//p_plus
                q_minus = N//p_minus

                if p_plus != 1 and p_plus != N:
                    p,q = p_plus, q_plus
                    break
                elif p_minus != 1 and p_minus != N:
                    p,q = p_minus, q_minus
                    break

            log['measure'].append(log_m0)
                    
            if p is not None and q is not None:
                r['factor'] = [p if p < q else q, q if q > p else p]
                break
        
        r['log'].append(log)

    print('---------------------------')
    print('answer:', N, '=', p if p < q else q, 'x', q if q > p else p)
    
    return r

factor_res = []
for i in range(1):
    print('\nRound '+str(i+1)+'\n')
    res = factor(25)
    print(json.dumps(res, indent=2, sort_keys=False))
    factor_res.append(res)

print(json.dumps(factor_res, indent=2, sort_keys=False))