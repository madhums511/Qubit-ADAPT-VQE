#!/usr/bin/env python3

import numpy as np
import math
from scipy.io import FortranFile
import scipy.optimize as opt
import os
import subprocess
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Operator
from qiskit.primitives import Estimator
from qiskit_ibm_runtime import Options, Session, QiskitRuntimeService  #,Estimator 
import matplotlib.pyplot as plt
from qiskit.circuit import Parameter
from qiskit.algorithms.optimizers import Optimizer
import time


#FUNCTIONS#############################

#runs fortran code and returns H matrix for given input
def h_eigvals(inp_):
    '''
    Params:
        inp_ : Input (formatted as per manual) to run the Fortran code to generate the Hamiltonian 
    Returns:
        QHQ_mat : H matrix of size corresponding to the input basis
    '''
    inp_list = inp_.split(' ')
    inp_list[3] = 'R' + inp_list[3]
    inp_list[4] = '_' + inp_list[4]

    #os.system('../H-H/hh/hhsp02.csh {inpu}'.format(inpu = inp_))
    subprocess.call(['hhsp02.csh', 'hh', 'b501', 'v08c', '7', '00'], stdout = subprocess.DEVNULL , stderr = subprocess.STDOUT)
    fname = 'scr/{filename}.hmt'.format(filename = (''.join(inp_list)))
    f = FortranFile(fname)

    #path = f'../H-H/hh/scr/hhb500vgo2R2_80.hmt'
    #f = FortranFile(path)

    h_triangular = f.read_reals(float)

    #print(h_triangular)

    #bsize = int(0.5 + math.sqrt(0.25 + 2 * len(h_triangular))) - 1 # the original basis size

    n = int(inp_list[1][1::])
    cdt = np.dtype(np.cdouble)
    QHQ_mat = np.zeros((n, n), dtype = cdt) # the Hamiltonian matrix

    #print(len(h_triangular))

    for j in np.arange(0, n):       #following fortran column-first indexing
        for i in np.arange(0, j+1):
            placeholder = int(i + j*(j+1)/2)
            QHQ_mat[i, j] = h_triangular[placeholder]
            QHQ_mat[j, i] = QHQ_mat[i, j]

    #evals, evecs = np.linalg.eig(QHQ_mat)
    return QHQ_mat #, evecs, np.sort(np.real(evals))

#truncates H matrix based on number of qubits available
def h_mat_trunc(h_mat, Nq):
    '''
    Params:
        h_mat : (list) Complete H matrix as generated by h_eigvals(inp_)
        Nq : (int) Number of available qubits
    Returns:
        h_mat_trunc : (list) Truncated H matrix of dimensions (Nq, Nq)
    '''
    n = 2**Nq # Hilbert space dimension = 2^Nq
    h_mat_trunc = []
    for i in range(n):
        h_mat_trunc.append(h_mat[i][0:n])
    
    return h_mat_trunc #matrix truncated to (n,n)

#decomposes H matrix into matrix elements in terms of Pauli strings
def gate_construction(Nq, row, column):
    '''
    Params:
        Nq : (int) Number of available qubits
        row : (int) Row index of the truncated H matrix
        column : (int) Column index of the truncated H matrix
    Returns:
        circ_list : (list(QuantumCircuit)) List of circuits of qubitised Hamiltonian framework
        sign_list : (list(int)) List of 0s and 1s indicating the sign of the corresponding expectation value of the circuit in circ_list (0 : + ; 1 : -)
    '''
    n = 2**Nq

    row_num = f'{row:0{Nq}b}' 
    col_num = f'{column:0{Nq}b}' 

    row_num = str(row_num)[::-1]
    col_num = str(col_num)[::-1]

    circ_list = [] 
    sign_list = np.zeros(n , dtype = int) 

    for i in range(n):
        circuit_init = QuantumCircuit(Nq)
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
                    circ_list[p].id(i)
                    p+=1
                    
                elif p < end:
                    circ_list[p].z(i)
                    p+=1
                else:
                    start = end
                    end += 2*step
                    mid = start +step
        
        elif (row_num[i] , col_num[i]) == ('1','1'):
            p = 0
            while end <= n:
                if p < mid:
                    circ_list[p].id(i)
                    p+=1
                elif p < end:
                    circ_list[p].z(i)
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
                    circ_list[p].x(i)
                    p+=1
                elif p < end:
                    circ_list[p].ry(((-1)*np.pi), i)
                    p+=1
                else:
                    start = end
                    end += 2*step
                    mid = start +step
        
        else:
            p = 0
            while end <= n:
                if p < mid:
                    circ_list[p].x(i)
                    p+=1
                elif p < end:
                    circ_list[p].ry(((-1)*np.pi), i)
                    sign_list[p] = int(sign_list[p]^1)
                    p+=1
                else:
                    start = end
                    end += 2*step
                    mid = start +step

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

#initialises state as Hartree Fock state where num_orb = Nq
def hf_state(num_orb):
    '''
    Params:
        num_orb : (int) Number of orbitals in representative molecule (equal to the number of available qubits in the system)
    Returns:
        qc : (QuantumCircuit) Circuit representing initial HF state
    '''
    qc = QuantumCircuit(num_orb)
    qc.x(0)
    qc.x(math.floor(num_orb/2))
    qc.barrier()
    return qc

#implements pauli string in circuit
def pauli_string(Nq, pool_element):
    string = pool_element[::-1]
    qcp = QuantumCircuit(Nq)
    qcp.barrier()
    for i in range(len(string)):
        if string[i] == 'Z':
            qcp.z(i)
        elif string[i] == 'Y':
            qcp.ry((-1)*np.pi, i)
    qcp.barrier()

    return qcp

#create exponential operator from pool using cnot staircase approach
def exp_V_stair(theta, Nq, pool_element):
    '''
    Params:
        theta : (Parameter) Parameter for pool element
        pool_element : (string) Pauli string
    Returns:
        qc : (QuantumCircuit) Circuit representing exponentiated pauli string parameterised by theta
    '''
    qc = QuantumCircuit(Nq)
    #qc.barrier()
    def cnot_stair(start_ind, end_ind, qc):
        j = start_ind
        while j < end_ind-1:
            p = j+ 1
            while pool_element[p] == 'I':
                p+=1
            qc.cx(j, p)
            j = p
        return qc
    
    def cnot_stair_rev(start_ind, end_ind, qc):
        j = start_ind
        while j > end_ind:
            p = j - 1
            while pool_element[p] == 'I':
                p-=1
            qc.cx(p,j)   
            j = p
        return qc 

    y_ind = []
    rz_ind = pool_element.find('Z')

    if rz_ind == -1:
        rz_ind = pool_element.find('Y')

    pool_element = pool_element[::-1]
    start_ind = pool_element.find('Y')

    for j in range(Nq - rz_ind):
        if pool_element[j] == 'Y':
            y_ind.append(j)
            qc.rz((-1)*(np.pi)/2, j)
            qc.h(j)

    qc_temp = QuantumCircuit(Nq)
    qc_temp = cnot_stair(start_ind, Nq-rz_ind, qc_temp)

    qc = qc.compose(qc_temp)
    qc.rz((-2)*theta, Nq-rz_ind-1)

    qc_temp2 = QuantumCircuit(Nq)
    qc_temp2 = cnot_stair_rev(Nq-rz_ind-1, start_ind, qc_temp2)
    qc = qc.compose(qc_temp2)

    for j in y_ind:
        qc.h(j)
        qc.rz((np.pi)/2, j)
        
    #qc.barrier()

    return qc

#applies ansatz to state
def ansatz(theta_num, param_list, Nq, pool_element, state, backend):
    '''
    Params:
        theta_num : (int) Parameter index
        param_list : (list(Parameter)) List storing indexed parameters
        pool_element : (string) Pauli string
        state : (QuantumCircuit) State to which new ansatz operator needs to be applied 
        backend: Runtime server according to which transpilation is performed
    Returns:
        transpile(qc, backend) : (QuantumCircuit) Transpiled circuit representing state with new ansatz operator appended
        param_list : (list(Parameter)) Updated list of indexed parameters
    '''
    param_list.append(Parameter('theta{num}'.format(num = theta_num)))
    print("Selected pool element: " , pool_element)
    qc = state.compose(exp_V_stair(param_list[theta_num], Nq, pool_element))
    return transpile(qc, backend), param_list

#Evaluates expectation value of operator H with respect to state psi
def expec_eval(h, psi, sim, hack):
    '''
    Params:
        h : (QuantumCircuit) Term in qubitised hamiltonian expansion
        psi : (QuantumCircuit) State for which expectation value of H is to be calculated
        sim : Estimator class- evaluates expectation value exactly (if imported from qiskit.primitives) or shot-by-shot (if imported from qiskit.runtime)
        hack : (boolean) if True uses matrix multipication instead of estimator class
    Returns:
        result.values[0] : (float) Expectation value of operator H wrt to state psi
    '''
    if not hack:
        op = Operator.from_circuit(h)
        job = sim.run(psi, op)
        result = job.result()
        return result.values[0]
    else:
        op = Operator(h)
        h_mat = op.to_matrix()
        psi_op = Operator(psi)
        psi_mat = psi_op.to_matrix()
        mat = (np.conj(psi_mat).T) @ h_mat @ psi_mat
        return np.real(mat[0][0])


#Evaluates expectation value of qubitised Hamiltonian H
def H_expec(H, psi, Nq, sim, hack):
    '''
    Params:
        H : (list) Truncated hamiltonian generated by the Fortran code
        psi : (QuantumCircuit) State for which expectation value of the hamiltonian is to be computed
        sim : Estimator class- evaluates expectation value exactly (if imported from qiskit.primitives) or shot-by-shot (if imported from qiskit.runtime)
    Returns:
        np.sum(temp1) : (float) Expectation value of the hamiltonian for state psi
    '''
    n = 2**Nq
    temp1 = []
    for row in range(n):
        temp2 = []
        for col in range(n):
            circ_list, sign_list = gate_construction(Nq, row, col)
            temp3 = []
            for k in range(len(circ_list)):
                temp3.append((expec_eval(circ_list[k], psi, sim, hack))*((-1)**(sign_list[k])))
            temp2.append((np.sum(temp3)/n)*(H[row][col]))
        temp1.append(np.sum(temp2))

    return np.real(np.sum(temp1))

#Evaluate commutator [H, A_i] for A_i \in {V}_n
def H_gradient(H, pool_element, psi, Nq, sim, thetas, hack): 
    '''
    Params:
        H : (list) Truncated hamiltonian generated by the Fortran code
        pool_element : (QuantumCircuit) Circuit representing a pool operator
        state : State for which expectation value of the commutator is to be computed
        thetas : (dict) Dictionary mapping circuit parameters to optimized values
    Returns:
        np.sum(temp1a)-np.sum(temp1b) : (float) Expectation value of the commutator for state 
    '''
    n = 2**Nq
    temp1 = []
    for row in range(n):
        temp2 = []
        for col in range(n):
            circ_list, sign_list = gate_construction(Nq, row, col)
            temp3 = []
            for k in range(len(circ_list)):
                temp3.append((expec_eval(circ_list[k].compose(pool_element), psi, sim, hack))*((-1)**(sign_list[k])))
            temp2.append((np.sum(temp3)/n)*(H[row][col]))
        temp1.append(np.sum(temp2))
    return -2*np.real(np.sum(temp1))

#Picks pool element with highest gradient and appends it to ansatz with appropriate variational parameter after checking for convergence
def ansatz_update(run_num, H, theta_list, pool_list, state, Nq, sim, conv_lim, param_list, V_list, hack):
    '''
    Params:
        run_num : (int) Number of iterations of the algorithm so far (equal to number of pool operators added to ansatz)
        H : (list) Truncated hamiltonian generated by the Fortran code
        theta_list (list(float)) List of optimised parameters
        pool_list : (list(string)) List of pauli strings for the pool of choice
        state : (QuantumCircuit) State for which a new ansatz operator needs to be applied
        conv_lim : (float) Convergence limit for gradient of expectation value
        param_list : (list(Parameter)) List of indexed parameters
    Returns:
        conv : (bool) True if convergence condition is met else False
        theta_list : (list(float)) Updated list of optimised parameters 
        transpile(qc, backend) : (QuantumCircuit) Transpiled circuit representing state with new ansatz operator appended
        param_list : (list(Parameter)) Updated list of indexed parameters
        V_list : (list(string)) List to keep track of appended pool elements
    '''
    grad_list = []
    pool_circ = []
    theta_init = 0 

    thetas = {param_list[i]: theta_list[i] for i in range(len(theta_list))}
    for p_ind in range(len(pool_list)):
        pool_circ.append(pauli_string(Nq, pool_list[p_ind]))
        psi = state.bind_parameters(thetas)
        grad_list.append(H_gradient(H, pool_circ[p_ind], psi, Nq, sim, thetas, hack)) #initialisation matters when calculating gradient (?)
    
    print("Gradients:" , {pool_list[i]:grad_list[i] for i in range(len(pool_list))})
    
    if np.linalg.norm(grad_list)<=conv_lim and run_num>1:
        return True, theta_list, V_list, (state, param_list)
    else:
        theta_list.append(theta_init)
        grad_list = list(np.abs(grad_list))
        element = pool_list[grad_list.index(max(grad_list))]
        V_list.append(element)
        return False, theta_list, V_list, ansatz(len(theta_list)-1, param_list, Nq, element , state, backend)
    

def minim_func(theta_list, param_list, psi, Nq, hack, sim):
    thetas = {param_list[i]: theta_list[i] for i in range(len(theta_list))}
    psi_ = psi.bind_parameters(thetas)
    energy = H_expec(H, psi_, Nq, sim, hack)
    return np.real(energy)

def grad_func(theta_list, param_list, psi, Nq, hack, sim, V_list):
    grad_list = np.zeros(len(V_list))
    pool_circ = []
    thetas = {param_list[i]: theta_list[i] for i in range(len(theta_list))}
    psi_ = psi.bind_parameters(thetas)
    pool_circ.append(pauli_string(Nq, V_list[-1]))
    grad_list[-1] = H_gradient(H, pool_circ[-1], psi_, Nq, sim, thetas, hack)
    
    for i in range(1, len(V_list)):
        rem = QuantumCircuit(Nq)
        for j in range(i):
            rem = rem.compose(exp_V_stair(param_list[-(j+1)], Nq, V_list[-(j+1)]))
        theta_rem1 = {param_list[-(k+1)] : -theta_list[-(k+1)] for k in range(i)}
        theta_rem2 = {param_list[-(k+1)] : theta_list[-(k+1)] for k in range(i)}
        op = (rem.bind_parameters(theta_rem1)).compose(pauli_string(Nq, V_list[-(i+1)]))
        op = op.compose(rem.bind_parameters(theta_rem2))
        grad_list[-(i+1)] = H_gradient(H, op, psi_, Nq, sim, thetas, hack)
        
    return list(np.real(grad_list))
#MAIN##################################

#######INPUTS##########################

inp = 'hh b501 v08c 7 00' 
Nq = 3 #Number of qubits in register
e_i = 0 #Eigenvector or eigenvalue index (for test1)
n = 2**Nq
hack = True #If hack is True, uses explicit matrix multiplication in place of the Estimator class
conv_lim = 1.e-06
conv = False
run_lim = 100

service= QiskitRuntimeService()
backend_name = 'ibmq_qasm_simulator'
backend = service.backend(backend_name)
options = Options()
options.resilience_level = 0
options.optimization_level = 0
options.execution.shots = 5
sim = Estimator() #when using qiskit.primitives


#Step 1: Qubitise Hamiltonian
H = h_mat_trunc(h_eigvals(inp), Nq)


#Step 2: Generate pool (pool1: Global pool, pool2: local pool)
pool_list = pool1(Nq)

#Step 3: Initialize state as HF state
state = QuantumCircuit(Nq) #hf_state(Nq)

theta_list = []
param_list = []
V_list = []
start = time.time()

energy_list = []
run_ind = 0

print("-----------------------------")
while run_ind < run_lim and not conv: 
    print('\n' , f'RUN {run_ind}' , '\n' )
    minim_ops = {"maxiter" : 300 , "disp" : True}
    conv, theta_list, V_list, (state, param_list) = ansatz_update(run_ind, H, theta_list, pool_list, state, Nq, sim, conv_lim, param_list, V_list, hack)
    result = opt.minimize(lambda x : minim_func(x, param_list, state, Nq, hack, sim), theta_list, method='BFGS', jac = lambda x : grad_func(x, param_list, state, Nq, hack, sim, V_list), tol = 1e-7,  options = minim_ops)
    theta_list = result.x
    theta_list = list(theta_list)
    print(f'Optimised parameters: {theta_list}')
    print("-----------------------------")
    run_ind+=1

thetas = {param_list[i]: theta_list[i] for i in range(len(theta_list))}
psi = state.bind_parameters(thetas)
energy = H_expec(H, psi, Nq, sim, hack)

evals, evecs = np.linalg.eig(H)
e0 = np.real(min(evals))
print("Lowest energy eigenvalue: ", e0)
print("Computed energy: ", np.real(energy))
print(f'Accuracy: {np.abs((e0 - np.real(energy))/e0)*100}%')
print("Time taken: " , time.time()-start)

plt.title(f'HaH: {inp} with {Nq} Qubits')
plt.plot(energy_list)
plt.plot(np.ones(len(energy_list))*e0)
plt.legend(["VQE energy" , f'Actual eigenvalue: {e0}'])