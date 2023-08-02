#!/usr/bin/env python3

import numpy as np
import math
from scipy.io import FortranFile
from scipy.optimize import minimize
from scipy.linalg import norm
import os
import subprocess
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Operator, Pauli
import matplotlib.pyplot as plt
from qiskit.algorithms.optimizers import Optimizer
from decimal import *
import time


#FUNCTIONS#############################

#runs fortran code and returns H matrix for given input
def h_mat(inp_ , Nq = None):
    '''
    Params:
        inp_ : Input (formatted as per manual) to run the Fortran code to generate the Hamiltonian 
    Returns:
        h_mat_trunc : H matrix of size corresponding to the input basis
    '''
    inp_list = inp_.split(' ')
    inp_list[3] = 'R' + inp_list[3]
    inp_list[4] = '_' + inp_list[4]

    #os.system('../H-H/hh/hhsp02.csh {inpu}'.format(inpu = inp_))
    subprocess.call(['hhsp02.csh', 'hh', 'b501', 'v08c', '7', '00'], stdout = subprocess.DEVNULL, stderr= subprocess.STDOUT)
    fname = 'scr/{filename}.hmt'.format(filename = (''.join(inp_list)))
    f = FortranFile(fname)

    #path = f'../H-H/hh/scr/hhb500vgo2R2_80.hmt'
    #f = FortranFile(path)

    h_triangular = f.read_reals(float)

    #print(h_triangular)

    #bsize = int(0.5 + math.sqrt(0.25 + 2 * len(h_triangular))) - 1 # the original basis size

    n = int(inp_list[1][1::])
    cdt = np.dtype(np.cdouble)
    QHQ_mat = np.zeros((n, n), dtype = float) # the Hamiltonian matrix

    #print(len(h_triangular))

    for j in np.arange(0, n):       #following fortran column-first indexing
        for i in np.arange(0, j+1):
            placeholder = int(i + j*(j+1)/2)
            QHQ_mat[i, j] = h_triangular[placeholder]
            QHQ_mat[j, i] = QHQ_mat[i, j]

    if Nq != None:
        n = 2**Nq # Hilbert space dimension = 2^Nq
        h_mat_trunc = []
        for i in range(n):
            h_mat_trunc.append(QHQ_mat[i][0:n])
    
        return h_mat_trunc 
    
    else:
        return QHQ_mat


#decomposes H matrix into matrix elements in terms of Pauli strings
def gate_construction(Nq, row, column):
    n = 2**Nq

    row_num = f'{row:0{Nq}b}' 
    col_num = f'{column:0{Nq}b}' 

    row_num = str(row_num)[::-1]
    col_num = str(col_num)[::-1]

    circ_list = [] 
    sign_list = np.zeros(n , dtype = int) 

    for i in range(n):
        circuit_init = ''
        circ_list.append(circuit_init)

    for i in range(Nq):
        end = n/(2**i)
        start = 0
        mid = end/2
        step = end/2
        #print('\n' ,row_num[i] , col_num[i])
        if (row_num[i] , col_num[i]) == ('0','0'):
            p = 0
            while end <= n:
                if p < mid:
                    circ_list[p] = 'I' + circ_list[p]
                    p+=1
                    
                elif p < end:
                    circ_list[p] = 'Z' + circ_list[p]
                    p+=1
                else:
                    start = end
                    end += 2*step
                    mid = start + step
        
        elif (row_num[i] , col_num[i]) == ('1','1'):
            p = 0
            while end <= n:
                if p < mid:
                    circ_list[p] = 'I' +circ_list[p]
                    p+=1
                elif p < end:
                    circ_list[p] = 'Z' +circ_list[p]
                    sign_list[p] = int(sign_list[p]^1)
                    p+=1
                else:
                    start = end
                    end += 2*step
                    mid = start +step

        elif (row_num[i] , col_num[i]) == ('0','1'):
            p = 0
            while end <= n:
                if p < mid:
                    circ_list[p] = 'X' +circ_list[p]
                    p+=1
                elif p < end:
                    circ_list[p] = 'Y' + circ_list[p]
                    p+=1
                else:
                    start = end
                    end += 2*step
                    mid = start +step
        
        else:
            p = 0
            while end <= n:
                if p < mid:
                    circ_list[p] = 'X' +circ_list[p]
                    p+=1
                elif p < end:
                    circ_list[p] = 'Y' +circ_list[p]
                    sign_list[p] = int(sign_list[p]^1)
                    p+=1
                else:
                    start = end
                    end += 2*step
                    mid = start +step
        
    for i in range(len(circ_list)):
        c = circ_list[i].count('Y')
        if c % 4 == 1:
            circ_list[i] = 'i' + circ_list[i]
        elif c % 4 == 2:
            circ_list[i] = '-' + circ_list[i]
        elif c % 4 == 3:
            circ_list[i] = '-i' + circ_list[i]

    return circ_list, sign_list

#global pool list
def pool1(Nq):
    '''
    Params:
        Nq : (int) Number of available qubits
    Returns:
        pool_list : (list(QuantumCircuit)) Global pool {V}_{Nq} list with 2*(Nq-1) elements
    '''
    pool_list = ['Y']
    for i in range(1 , Nq):
        for j in range(len(pool_list)):
            pool_list[j]+='Z'
        pool_list.append('Y')

        if i>1:
            pool_list.append("YI")
        
    for i in range(len(pool_list)):
        pool_list[i] = pool_list[i][::-1]+ ('I'*(Nq - len(pool_list[i]))) 
    return pool_list

#local pool list
def pool2(Nq):
    '''
    Params:
        Nq : (int) Number of available qubits
    Returns:
        pool_list : (list(QuantumCircuit)) Local pool {G}_{Nq} list with 2*(Nq-1) elements
    '''
    pool_list = []
    for i in range(Nq-1):
        pool_list.append('I'*i)
        pool_list[i]+='ZY'
        pool_list[i]+= ('I'*(Nq-2-i))
    for j in range((Nq-1), 2*(Nq-1)):
        s = j - (Nq-1)
        pool_list.append('I'*s)
        pool_list[j]+='Y'
        pool_list[j]+= ('I'*(Nq-1-s))
    return pool_list


#implements pauli string in circuit
def pauli_string(pool_element):
    pauli_op = Pauli(pool_element)
    pauli_mat = pauli_op.to_matrix()
    return pauli_mat

#create exponential operator from pool using cnot staircase approach
def exp_mat(theta, Nq, pool_element):
    '''
    Params:
        theta : (Parameter) Parameter for pool element
        pool_element : (string) Pauli string
    Returns:
        qc : (QuantumCircuit) Circuit representing exponentiated pauli string parameterised by theta
    '''
    pauli_mat = pauli_string(pool_element)
    pauli_id = np.identity(2**Nq)
    return (np.cos(theta))*pauli_id + (1.0j)*(np.sin(theta))*pauli_mat

#applies ansatz to state
def ansatz(theta_list, Nq, V_list):
    state = np.identity(2**Nq)
    for i in range(len(V_list)):
        state = exp_mat(theta_list[i] , Nq, V_list[i]) @ state

    return state

#Evaluates expectation value of operator H with respect to state psi
def H_expec(H, psi, Nq):
    n = 2**Nq
    matsum = 0
    for row in range(n):
        for col in range(n):
            circ_list, sign_list = gate_construction(Nq, row, col)
            for i in range(len(circ_list)):
                bra = H[row][col]*(pauli_string(circ_list[i]) @ psi)
                mat = ((np.conj(bra).T) @  psi)
                matsum += mat[0][0]*((-1)**sign_list[i])
    return np.real(matsum)/n

def H_gradient(H, psi, pool_element, Nq):
    n = 2**Nq
    matsum = 0
    for row in range(n):
        for col in range(n):
            circ_list, sign_list = gate_construction(Nq, row, col)
            for i in range(len(circ_list)):
                bra = H[row][col]*(pauli_string(circ_list[i]) @ psi)
                ket = (pauli_string('i'+ pool_element) @  psi)
                mat = (np.conj(bra).T) @ ket
                matsum += mat[0][0]*((-1)**sign_list[i])

    return 2*(np.real(matsum))/n


#Picks pool element with highest gradient and appends it to ansatz with appropriate variational parameter after checking for convergence
def V_update(run_num, H, theta_list, pool_list, Nq, conv_lim, V_list):
    grad_list = []
    theta_init = 0.0
    psi = ansatz(theta_list, Nq, V_list)
    for p_ind in range(len(pool_list)):
        grad_list.append(H_gradient(H, psi, pool_list[p_ind], Nq)) #initialisation matters when calculating gradient (?)
    
    print("Gradients:" , {pool_list[i]:grad_list[i] for i in range(len(pool_list))})
    
    if norm(grad_list)<=conv_lim and run_num>1:
        return True, theta_list, V_list
    else:
        theta_list.append(theta_init)
        grad_list = list(np.abs(grad_list))
        element = pool_list[grad_list.index(max(grad_list))]
        V_list.append(element)
        print(f'Selected operators: {V_list}')
        return False, theta_list, V_list
    

#Minimizes theta where psi is the state with the updated ansatz (using bfgs)
    
def minim_func(theta_list , Nq, V_list, H):
    psi = ansatz(theta_list, Nq, V_list)
    energy = H_expec(H, psi, Nq)
    #energy_list.append(energy)
    return energy

'''
def grad_func(theta_list , Nq, V_list, H):
    grad_list = []
    psi = ansatz(theta_list, Nq, V_list)
    for i in range(len(V_list)):
        grad_list.append(H_gradient(H, psi, V_list[i] ,Nq))
    print(f"Grad val: {theta_list}, { list(np.real(grad_list))}")
    return list(np.real(grad_list))
'''

def grad_func(theta_list, Nq, V_list, H):
    grad_list = np.zeros(len(V_list))
    psi = ansatz(theta_list, Nq, V_list)
    grad_list[-1] = H_gradient(H, psi, V_list[-1], Nq)
    bra = (H @ ansatz(theta_list, Nq, V_list))
    for i in range(1, len(V_list)):
        rem = np.conj(exp_mat(theta_list[-i], Nq, V_list[-i])).T
        bra = rem @ bra
        psi = rem @ psi
        mat = np.conj(bra).T  @ pauli_string('i'+ V_list[-(i+1)]) @ psi
        grad_list[-(i+1)] = 2*np.real(mat[0][0])
    #print(f'Grad val: {theta_list} {np.real(grad_list)}')
    return list(np.real(grad_list))


#MAIN##################################

#######INPUTS##########################

inp = 'hh b501 v08c 7 00' 
Nq = 4 #Number of qubits in register
e_i = 0 #Eigenvector or eigenvalue index (for test1)
n = 2**Nq
hack = True #If hack is True, uses explicit matrix multiplication in place of the Estimator class
conv_lim = 1.e-06
conv = False
run_lim = 2**Nq
minim_ops = {"maxiter" : 300 , "disp" : True }

#Step 1: Qubitise Hamiltonian
H = h_mat(inp, Nq = Nq)

#Step 2: Generate pool (pool1: Global pool, pool2: local pool)
pool_list = pool1(Nq)

#Step 3: Initialize 
theta_list = []
param_list = []
V_list = []
start = time.time()

run_ind = 0

print("-----------------------------")
while run_ind < run_lim and not conv: 
    print('\n' , f'RUN {run_ind}' , '\n' )
    conv, theta_list, V_list = V_update(run_ind, H, theta_list, pool_list, Nq, conv_lim, V_list)
    result = minimize(lambda x: minim_func(x , Nq, V_list, H), theta_list, method='BFGS', jac= lambda x: grad_func(x , Nq, V_list, H), tol = 1.e-7, options = minim_ops)
    theta_list = list(result.x)
    print(f'Optimised parameters: {theta_list}')
    print("-----------------------------")
    run_ind+=1

f_energy = minim_func(theta_list, Nq, V_list, H)
evals, evecs = np.linalg.eig(H)
e0 = np.real(min(evals))
print("Lowest energy eigenvalue: ", e0)
print("Computed energy: ", np.real(f_energy))
print(f'Accuracy: {np.abs((e0 - np.real(f_energy))/e0)*100}%')
print("Time taken: " , time.time()-start)

'''
plt.title(f'HaH: {inp} with {Nq} Qubits')
plt.plot(f_energy)
plt.plot(np.ones(len(energy_list))*e0)
plt.legend(["VQE energy" , f'Actual eigenvalue: {e0}'])
plt.show()
'''
exit()