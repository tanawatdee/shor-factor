#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'inline')

import math
import random
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, execute, Aer, IBMQ
from qiskit.tools.jupyter import *
from qiskit.visualization import *

provider = IBMQ.load_account()
backend = Aer.get_backend('qasm_simulator')

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

class CR(ClassicalRegister):
    def __init__(self, *args, **kwargs):
        super(CR, self).__init__(*args, **kwargs)

class QR(QuantumRegister):
    def __init__(self, *args, **kwargs):
        super(QR, self).__init__(*args, **kwargs)

class QC(QuantumCircuit):
    def __init__(self, *args, **kwargs):
        super(QC, self).__init__(*args, **kwargs)
        
    def plot_prob(self):
        job = execute(self, backend, shots=0x100000)
        job.result().get_counts(qc)
        return plot_histogram(job.result().get_counts(qc), figsize=[50,10])
    
    def get_measure(self, shots):
        job = execute(self, backend, shots=shots, memory=True)
        return job.result().get_memory()
    
    ###################################################
    # constant - set qubits from |0> to |c>
    # c - a constant to set
    # q - qubits to set           n qubits
    # ctrl - control (optional)   1 qubit
    ###################################################
    def constant(self, c, q, ctrl=None):
        self.barrier()
        for i in range(len(q)):
            self.x(q[i]) if c>>i&1 and ctrl is None else self.cx(ctrl, q[i]) if c>>i&1 else 0
        self.barrier()
            
    ###################################################
    # add - Controlled adder
    #              adding two numbers a + b = s
    # ctrl - control    1 qubit
    # a -> a            n qubits
    # b -> s            n qubits
    # z -> z xor s_n+1  1 qubit
    ###################################################
    def add(self, ctrl, a, b, z):
        self.ccx(ctrl, a[-1], z)
        for i in range(1, len(a)):
            self.cx(a[-i], b[-i])
            self.barrier()
        for i in range(2, len(a)):
            self.cx(a[-i], a[-i+1])
        for i in range(0, len(a)-1):
            self.ccx(a[i], b[i], a[i+1])
        for i in range(2):
            self.ccx(a[-1], b[-1], b[0])
            self.ccx(ctrl, b[0], z)
        for i in range(1, len(a)):
            self.ccx(ctrl, a[-i], b[-i])
            self.ccx(a[-i-1], b[-i-1], a[-i])
        for i in range(1, len(a)-1):
            self.cx(a[i], a[i+1])
            self.barrier()
        for i in range(1, len(a)):
            self.cx(a[i],  b[i])
        self.ccx(ctrl, a[0], b[0])
        self.barrier()
        
    ###################################################
    # const_add - Controlled constant adder
    #                  adding two numbers a + const
    # const - a constant to be added
    # ctrl - controls   1 qubit or list of qubit
    # a -> a + const    n qubits
    ###################################################
    def const_add(self, const, ctrl, a):
        ctrl = ctrl if type(ctrl) == list else [ctrl]
        const_bin = list(map(lambda i: const>>i&1, range(math.ceil(math.log(abs(const)+1, 2)) + 1)))
        
        for i in range(2, len(a)+1):
            for j in range(1, i):
                self.mcrx(math.pi/2**(i-j), [a[-i]], a[-j])
            self.barrier()
        for i in range(len(a)):
            theta = math.pi*sum(map(lambda j: const_bin[j]/2**(i-j), range(i+1)))
            self.mcrx(theta, ctrl, a[i])
        self.barrier()
        for i in range(len(a), 1, -1):
            for j in range(i-1, 0, -1):
                self.mcrx(-math.pi/2**(i-j), [a[-i]], a[-j])
            self.barrier()
            
    ###################################################
    # cmp - Controlled comparator
    #            flip bit z if a <= b
    # ctrl - control       1 qubit
    # a -> a               n qubits
    # b -> b               n qubits
    # z -> z xor (a <= b)  1 qubit
    ###################################################
    def cmp(self, ctrl, a, b, z):
        for i in range(len(a)-1, -1, -1):
            self.cx(a[i],  b[i])
        for i in range(len(a)-2, 0, -1):
            self.cx(a[i], a[i+1])
            self.barrier()
        for i in range(len(a)-1, 0, -1):
            self.ccx(a[-i-1], b[-i-1], a[-i])
            self.cx(a[-i], b[-i])
        for i in range(2):
            self.ccx(a[-1], b[-1], b[0])
            self.ccx(ctrl, b[0], z)
        for i in range(len(a)-2, -1, -1):
            self.ccx(a[i], b[i], a[i+1])
        for i in range(len(a)-1, 1 ,-1):
            self.cx(a[-i], a[-i+1])
        
        self.ccx(ctrl, a[-1], z)
        # reverse
        for i in range(2, len(a)):
            self.cx(a[-i], a[-i+1])
        for i in range(0, len(a)-1):
            self.ccx(a[i], b[i], a[i+1])
        for i in range(1, len(a)):
            self.cx(a[-i], b[-i])
            self.ccx(a[-i-1], b[-i-1], a[-i])
        for i in range(1, len(a)-1):
            self.cx(a[i], a[i+1])
            self.barrier()
        for i in range(len(a)):
            self.cx(a[i],  b[i])
        self.cx(ctrl, z)
        self.barrier()
        
    ###################################################
    # add_mod - Controlled adder modulo m
    #              adding two numbers a + b = s
    # m - modulo
    # ctrl - control    1 qubit
    # a -> a            n qubits
    # b -> s            n qubits
    # z -> z xor s_n+1  1 qubit
    ###################################################
    def add_mod(self, m, ctrl, a, b, z):
        self.add(ctrl, a, b, z)
        self.const_add(-m, ctrl, [*b, z])
        self.const_add(m, [ctrl, z], b)
        self.cmp(ctrl, a, b, z)
    
    ###################################################
    # dbl - Controlled double using cyclic shift
    # ctrl - control    1 qubit
    # a -> 2*a          n qubits
    # z -> 2*a_n+1      1 qubit
    ###################################################
    def dbl(self, ctrl, a, z):
        self.cswap(ctrl, z, a[-1])
        for i in range(1, len(a)):
            self.cswap(ctrl, a[-i], a[-i-1])
            
    ###################################################
    # dbl_mod - Controlled double modulo m
    # m - modulo
    # ctrl - control    1 qubit
    # a -> 2*a          n qubits
    # z -> 2*a_n+1      1 qubit
    ###################################################
    def dbl_mod(self, m, ctrl, a, z):
        self.dbl(ctrl, a, z)
        self.const_add(-m, ctrl, [*a, z])
        self.const_add(m, [ctrl, z], a)
        self.x(a[0])
        self.ccx(ctrl, a[0], z)
        self.x(a[0])
        
    ###################################################
    # swap_group - Controlled swap between
    #                   2 groups of qubits
    # ctrl - control    1 qubit
    # a - 1st group     n qubits
    # b - 2nd group     n qubits
    ###################################################
    def swap_group(self, ctrl, a, b):
        for i in range(len(a)):
            self.cswap(ctrl, a[i], b[i])
            
    ###################################################
    # mul_mod - Controlled multiplication modulo m
    # m - modulo
    # const - constant to multiply
    # ctrl - control    1 qubit
    # a -> a            n qubits
    # b -> const*a%m    n qubits
    # z -> |0> ancilla  1 qubit
    ###################################################
    def mul_mod(self, m, const, ctrl, a, b, z):
        self.add_mod(m, ctrl, a, b, z)
        for i in range(math.ceil(math.log(const+1, 2)) -2, -1, -1):
            self.dbl_mod(m, ctrl, b, z)
            if const>>i&1:
                self.add_mod(m, ctrl, a, b, z)
    
    ###################################################
    # mul_mod_in - Controlled multiplication
    #                   modulo m (in-place)
    # m - modulo
    # const - constant to multiply
    # ctrl - control                   1 qubit
    # q_sup - control qubit regiter  2*n qubits
    # q_mod -> const*q_mod%m           n qubits
    # q_acc -> |0> ancillae            n qubits
    # q_anc -> |0> ancilla             1 qubit
    ###################################################
    def mul_mod_in(self, m, const, ctrl, q_sup, q_mod, q_acc, q_anc):
        const_inv = mod_inverse(m, const)
        qc_inv    = QC(q_sup, q_mod, q_acc, q_anc)
        qc_inv.mul_mod(m, const_inv, ctrl, q_mod, q_acc, q_anc[0])
        qc_inv = qc_inv.inverse()

        self.mul_mod(m, const, ctrl, q_mod, q_acc, q_anc[0])
        self.swap_group(ctrl, q_mod, q_acc)
        self.extend(qc_inv)
        
    ###################################################
    # EXP - Controlled modular exponentiation
    # m - modulo
    # const - constant to multiply
    # q_sup - control qubit regiter  2*n qubits
    # q_mod -> const**q_sup%m          n qubits
    # q_acc -> |0> ancillae            n qubits
    # q_anc -> |0> ancilla             1 qubit
    ###################################################
    def EXP(self, m, const, q_sup, q_mod, q_acc, q_anc):
        for i in range(len(q_sup)):
            self.mul_mod_in(m, const, q_sup[i], q_sup, q_mod, q_acc, q_anc)
            const = const**2%m
    
    ###################################################
    # QFT - Quantum Fourier Transform
    # q - input and output
    ###################################################
    def QFT(self, q):
        n = len(q)
        
        self.barrier()
        for i in range(n//2):
            self.swap(q[i], q[n-1-i])

        self.barrier()
        for i in range(n):
            for j in range(i):
                self.cu1(math.pi/2**(i-j), q[i], q[j])
                self.barrier()
            self.h(q[i])
            self.barrier()


# In[9]:


N = 21 # number to factor
N_BIT = math.ceil(math.log(N+1, 2)) # N is represented in N_BIT bits
N_QFT = 2**(2*N_BIT) # number of basis states in QFT
p,q = None,None # 1st factor, 2nd factor

while p is None or q is None:
    
    if N%2 == 0: # quantum part is not designed for even number
        p,q = 2, N//2
        print(N, 'is divisible by 2.')
        break
        
    X = random.randint(2, N-1) # chooe a random number 1<X<N
    
    print('---------------------------')
    print('choose X =', X)
    
    gcdXN = gcd(X,N)

    if gcdXN != 1:
        p,q = gcdXN, N//gcdXN
        print(X, 'and', N, 'are not coprime.')
        print('apply gcd(', X, ',', N, ')')
        break
    
    ############## This is quantum part #############
    
    q_sup = QR(2*N_BIT)  # uniform superposition qubits
    q_mod = QR(N_BIT)    # modular exponentiation qubits
    q_acc = QR(N_BIT)    # ancillae for accumulation
    q_anc = QR(1)        # ancillae
    c_qft = CR(2*N_BIT)  # QFT measurement result
    qc    = QC(q_sup, q_mod, q_acc, q_anc, c_qft)

    qc.h(q_sup)
    qc.x(q_mod[0])

    qc.EXP(N, X, q_sup, q_mod, q_acc, q_anc)
    qc.QFT(q_sup)

    qc.measure(q_sup, c_qft)

    m = qc.get_measure(2*N_BIT)

    #################### end of quantum part ##################

    for m0 in m:
        m0  = int(m0,2)
        seq = cont_frac(m0, N_QFT)
        print('')
        print('measured value:', m0, '(of', N_QFT, 'basis states)')
        print('continued fraction:', seq)

        for j in range(1, len(seq)):
            k0,r0 = simp_frac(seq[:j+1])

            if r0 > N:
                break

            if abs(k0/r0 - m0/N_QFT) > 1/N_QFT:
                continue

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

        if p is not None and q is not None:
            break

print('---------------------------')
print('answer:', N, '=', p if p < q else q, 'x', q if q > p else p)

