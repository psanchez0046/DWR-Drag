import numpy as np
import glob
import sys
import time 
from scipy import sparse
# from scipy.signal import find_peaks
import scipy.sparse.linalg
from ypstruct import struct
    
def solve_NS_DWR_ll(Re1, Re2, NN, Y, H, G1, G2, R1, R2, R6, ringW, stepW, ringSubs, upperBC, R5_adim):
    # This function solves the Navier-Stokes equations with no-slip
    # and Boussinesq_Scriven boundary conditions with second order centered 
    # finite differences for the DWR configuration
    # INPUTS:
    # Re1 and Re2: Reynolds number
    # Y: "Boussinesq" number
    # NN: Bulk Viscosity ratio
    # Sres: mesh spacing
    # OUTPUT:
    # g: velocity field (nondimensional)
	
    Sres = ringW/ringSubs
    stepSubs = int(np.round(stepW/Sres))
    ringSubs_2 = int(ringSubs/2)

    NG1 = int(np.round(G1/Sres))
    NG2 = int(np.round(G2/Sres))
    N1 = NG1 + NG2 + ringSubs + 2*stepSubs
    N2 = NG1 + NG2 + ringSubs
    M = int(np.round(H/Sres))
    dim = (N1+1)*(M+1) + (N2+1)*M# Dimension of the square coefficients matrix A 
    # Adimensional parameters
    R2_adim = R2/R6
    R1_adim = R1/R6
    Sres_adim = Sres/R6
    AA = 1/Sres_adim

    ## ACUMULATING INDEX AND VALUES FOR A
    ## Upper boundary condition
    if upperBC == 'ns':
        # no-slip (closed contour)
        rows = np.array(np.arange(1, N1, dtype='int'))
        cols = np.array(np.arange(1, N1, dtype='int'))
        coefs = np.ones(N1-1, dtype='complex')
    elif upperBC == 'fb':
        # free interface nodes (open contour)
        rows = np.concatenate((np.arange(1, N1, dtype='int'),
                np.arange(1, N1, dtype='int'),
                np.arange(1, N1, dtype='int'),
                np.arange(1, N1, dtype='int')))
        cols = np.concatenate((np.arange(1, N1, dtype='int') - 1,
                np.arange(1, N1, dtype='int'),
                np.arange(1, N1, dtype='int') + 1,
                np.arange(1, N1, dtype='int') + (N1 + 1)))
        rind = np.arange(1, N1, dtype='int')
        coefs = np.concatenate((1 - 1/(2*(AA*R2_adim + rind)),
                    -1j*(Re2/(AA*AA)) - 4 - 1./((AA*R2_adim + rind)*(AA*R2_adim + rind)),
                    1 + 1./(2*(AA*R2_adim + rind)),
                    2*np.ones(N1-1, dtype='complex')))
    ## Walls
    # Upper left wall
    rows = np.append(rows, np.arange(M+1, dtype='int')*(N1+1))
    cols = np.append(cols, np.arange(M+1, dtype='int')*(N1+1))
    coefs = np.append(coefs, np.ones(M+1, dtype='complex'))
    # Upper right wall
    rows = np.append(rows, np.arange(1, M+2, dtype='int')*(N1+1)-1)
    cols = np.append(cols, np.arange(1, M+2, dtype='int')*(N1+1)-1)
    coefs = np.append(coefs, np.ones(M+1, dtype='complex'))
    # Lower left wall
    rows = np.append(rows, np.arange(M, dtype='int')*(N2+1)+(N1+1)*(M+1))
    cols = np.append(cols, np.arange(M, dtype='int')*(N2+1)+(N1+1)*(M+1))
    coefs = np.append(coefs, np.ones(M, dtype='complex'))
    # Lower right wall
    rows = np.append(rows, np.arange(1, M+1, dtype='int')*(N2+1)+(N1+1)*(M+1)-1)
    cols = np.append(cols, np.arange(1, M+1, dtype='int')*(N2+1)+(N1+1)*(M+1)-1)
    coefs = np.append(coefs, np.ones(M, dtype='complex'))
    # Ground
    rows = np.append(rows, np.arange(1, N2, dtype='int')+(M-1)*(N2+1)+(N1+1)*(M+1))
    cols = np.append(cols, np.arange(1, N2, dtype='int')+(M-1)*(N2+1)+(N1+1)*(M+1))
    coefs = np.append(coefs, np.ones(N2-1, dtype='complex'))
    # Left step
    rows = np.append(rows, np.arange(1, stepSubs+1, dtype='int')+M*(N1+1))
    cols = np.append(cols, np.arange(1, stepSubs+1, dtype='int')+M*(N1+1))
    coefs = np.append(coefs, np.ones(stepSubs, dtype='complex'))
    # Right step
    rows = np.append(rows, np.arange(stepSubs, dtype='int')+(M+1)*(N1+1)-stepSubs-1)
    cols = np.append(cols, np.arange(stepSubs, dtype='int')+(M+1)*(N1+1)-stepSubs-1)
    coefs = np.append(coefs, np.ones(stepSubs, dtype='complex'))
    ## Interface
    # Inner nodes
    rows = np.append(rows, [np.arange(1, NG1, dtype='int')+M*(N1+1)+stepSubs,
                            np.arange(1, NG1, dtype='int')+M*(N1+1)+stepSubs,
                            np.arange(1, NG1, dtype='int')+M*(N1+1)+stepSubs,
                            np.arange(1, NG1, dtype='int')+M*(N1+1)+stepSubs,
                            np.arange(1, NG1, dtype='int')+M*(N1+1)+stepSubs])
    cols = np.append(cols, [np.arange(1, NG1, dtype='int')+(M-1)*(N1+1)+stepSubs,
                            np.arange(1, NG1, dtype='int')+M*(N1+1)+stepSubs-1,
                            np.arange(1, NG1, dtype='int')+M*(N1+1)+stepSubs,
                            np.arange(1, NG1, dtype='int')+M*(N1+1)+stepSubs+1,
                            np.arange(1, NG1, dtype='int')+(M+1)*(N1+1)])
    rind = (np.arange(1, NG1, dtype='complex') + stepSubs)
    coefs = np.append(coefs, [(1/(AA*Y))*np.ones(NG1-1, dtype='complex'),
                                (1 - 1/(2*(AA*R2_adim + rind)))*(NN + (1/(2*AA))*(1 + 1/Y)),
                                -NN*(2 + 1/((AA*R2_adim + rind)*(AA*R2_adim + rind))) - (1/(2*AA))*(1j*(Re1/(AA*AA)) + 4 + 1/((AA*R2_adim + rind)*(AA*R2_adim + rind)) + (1/Y)*(1j*(Re2/(AA*AA)) + 4 + 1/((AA*R2_adim + rind)*(AA*R2_adim + rind)))),
                                (1 + 1./(2*(AA*R2_adim + rind)))*(NN + (1/(2*AA))*(1+1/Y)),
                                (1./(AA))*np.ones(NG1-1, dtype='complex')])
    # Outer nodes
    rows = np.append(rows, [np.arange(NG2-1, dtype='int')+NG1+M*(N1+1)+stepSubs+ringSubs+1,
                            np.arange(NG2-1, dtype='int')+NG1+M*(N1+1)+stepSubs+ringSubs+1,
                            np.arange(NG2-1, dtype='int')+NG1+M*(N1+1)+stepSubs+ringSubs+1,
                            np.arange(NG2-1, dtype='int')+NG1+M*(N1+1)+stepSubs+ringSubs+1,
                            np.arange(NG2-1, dtype='int')+NG1+M*(N1+1)+stepSubs+ringSubs+1])
    cols = np.append(cols, [np.arange(NG2-1, dtype='int')+NG1+(M-1)*(N1+1)+stepSubs+ringSubs+1,
                            np.arange(NG2-1, dtype='int')+NG1+M*(N1+1)+stepSubs+ringSubs,
                            np.arange(NG2-1, dtype='int')+NG1+M*(N1+1)+stepSubs+ringSubs+1,
                            np.arange(NG2-1, dtype='int')+NG1+M*(N1+1)+stepSubs+ringSubs+2,
                            np.arange(NG2-1, dtype='int')+NG1+(M+1)*(N1+1)+stepSubs+1])
    rind = np.arange(1, NG2, dtype='complex')+NG1+stepSubs+ringSubs
    coefs = np.append(coefs, [(1/(AA*Y))*np.ones(NG2-1, dtype='complex'),
                                (1 - 1/(2*(AA*R2_adim + rind)))*(NN + (1/(2*AA))*(1 + 1/Y)),
                                -NN*(2 + 1/((AA*R2_adim + rind)*(AA*R2_adim + rind))) - (1/(2*AA))*(1j*(Re1/(AA*AA)) + 4 + 1/((AA*R2_adim + rind)*(AA*R2_adim + rind)) + (1/Y)*(1j*(Re2/(AA*AA)) + 4 + 1/((AA*R2_adim + rind)*(AA*R2_adim + rind)))),
                                (1 + 1./(2*(AA*R2_adim + rind)))*(NN + (1/(2*AA))*(1+1/Y)),
                                (1./(AA))*np.ones(NG2-1, dtype='complex')])

    ## Upper Bulk
    P, Q = np.mgrid[1:M-ringSubs_2, 1:N1]
    rows = np.append(rows, [(P*(N1+1)+Q).reshape(1, -1),
                            (P*(N1+1)+Q).reshape(1, -1),
                            (P*(N1+1)+Q).reshape(1, -1),
                            (P*(N1+1)+Q).reshape(1, -1),
                            (P*(N1+1)+Q).reshape(1, -1)])
    cols = np.append(cols, [((P-1)*(N1+1)+Q).reshape(1, -1),
                            (P*(N1+1)+Q-1).reshape(1, -1),
                            (P*(N1+1)+Q).reshape(1, -1),
                            (P*(N1+1)+Q+1).reshape(1, -1),
                            ((P+1)*(N1+1)+Q).reshape(1, -1)])
    rind = np.arange(2, N1+1, dtype='complex')-1
    coefs = np.append(coefs, [np.ones((N1-1)*(M-ringSubs_2-1), dtype='complex'),
                                np.tile(1 - 1/(2*(AA*R2_adim + rind)), M-ringSubs_2-1),
                                np.tile(-1j*(Re2/(AA*AA)) - 4 - 1/((AA*R2_adim + rind)*(AA*R2_adim + rind)), M-ringSubs_2-1),
                                np.tile(1 + 1./(2*(AA*R2_adim + rind)), M-ringSubs_2-1),
                                np.ones((N1-1)*(M-ringSubs_2-1), dtype='complex')])
    for k in range(ringSubs_2):
        # Inner
        rows = np.append(rows, [np.arange(1, NG1+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k)*(N1+1),
                                np.arange(1, NG1+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k)*(N1+1),
                                np.arange(1, NG1+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k)*(N1+1),
                                np.arange(1, NG1+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k)*(N1+1),
                                np.arange(1, NG1+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k)*(N1+1)])
        cols = np.append(cols, [np.arange(1, NG1+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k-1)*(N1+1),
                                np.arange(1, NG1+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k)*(N1+1)-1,
                                np.arange(1, NG1+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k)*(N1+1),
                                np.arange(1, NG1+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k)*(N1+1)+1,
                                np.arange(1, NG1+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k+1)*(N1+1)])
        rind = np.arange(1, NG1+stepSubs+ringSubs_2-k, dtype='complex')
        coefs = np.append(coefs, [np.ones(NG1+stepSubs+ringSubs_2-k-1, dtype='complex'),
                                    1 - 1/(2*(AA*R2_adim + rind)),
                                    -1j*(Re2/(AA*AA)) - 4 - 1/((AA*R2_adim + rind)*(AA*R2_adim + rind)),
                                    1 + 1/(2*(AA*R2_adim + rind)),
                                    np.ones(NG1+stepSubs+ringSubs_2-k-1, dtype='complex')])
        # Outer
        rows = np.append(rows, [np.arange(1, NG2+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k)*(N1+1)+NG1+stepSubs+ringSubs_2+k,
                                np.arange(1, NG2+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k)*(N1+1)+NG1+stepSubs+ringSubs_2+k,
                                np.arange(1, NG2+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k)*(N1+1)+NG1+stepSubs+ringSubs_2+k,
                                np.arange(1, NG2+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k)*(N1+1)+NG1+stepSubs+ringSubs_2+k,
                                np.arange(1, NG2+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k)*(N1+1)+NG1+stepSubs+ringSubs_2+k])
        cols = np.append(cols, [np.arange(1, NG2+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k-1)*(N1+1)+NG1+stepSubs+ringSubs_2+k,
                                np.arange(1, NG2+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k)*(N1+1)+NG1+stepSubs+ringSubs_2+k-1,
                                np.arange(1, NG2+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k)*(N1+1)+NG1+stepSubs+ringSubs_2+k,
                                np.arange(1, NG2+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k)*(N1+1)+NG1+stepSubs+ringSubs_2+k+1,
                                np.arange(1, NG2+stepSubs+ringSubs_2-k, dtype='int')+(M-ringSubs_2+k+1)*(N1+1)+NG1+stepSubs+ringSubs_2+k])
        rind = np.arange(1, NG2+stepSubs+ringSubs_2-k, dtype='complex')+NG1+stepSubs+ringSubs_2
        coefs = np.append(coefs, [np.ones(NG2+stepSubs+ringSubs_2-k-1, dtype='complex'),
                                    1 - 1./(2*(AA*R2_adim + rind)),
                                    -1j*(Re2/(AA*AA)) - 4 - 1/((AA*R2_adim + rind)*(AA*R2_adim + rind)),
                                    1 + 1./(2*(AA*R2_adim + rind)),
                                    np.ones(NG2+stepSubs+ringSubs_2-k-1, dtype='complex')])
    # Lower Bulk
    for k in range(ringSubs_2):
        # Inner
        rows = np.append(rows, [np.arange(NG1+k, dtype='int')+(N1+1)*(M+1)+(N2+1)*k+1,
                                np.arange(NG1+k, dtype='int')+(N1+1)*(M+1)+(N2+1)*k+1,
                                np.arange(NG1+k, dtype='int')+(N1+1)*(M+1)+(N2+1)*k+1,
                                np.arange(NG1+k, dtype='int')+(N1+1)*(M+1)+(N2+1)*k+1,
                                np.arange(NG1+k, dtype='int')+(N1+1)*(M+1)+(N2+1)*k+1])
        cols = np.append(cols, [np.arange(NG1+k, dtype='int')+(N1+1)*M+stepSubs+(N2+1)*k+(k>0)*stepSubs+1,
                                np.arange(NG1+k, dtype='int')+(N1+1)*(M+1)+(N2+1)*k,
                                np.arange(NG1+k, dtype='int')+(N1+1)*(M+1)+(N2+1)*k+1,
                                np.arange(NG1+k, dtype='int')+(N1+1)*(M+1)+(N2+1)*k+2,
                                np.arange(NG1+k, dtype='int')+(N1+1)*(M+1)+(N2+1)*(k+1)+1])
        rind = np.arange(NG1+k)+1
        coefs = np.append(coefs, [np.ones(NG1+k, dtype='complex'),
                                    1 - 1/(2*(AA*R1_adim + rind)),
                                    -1j*(Re1/(AA*AA)) - 4 - 1/((AA*R1_adim + rind)*(AA*R1_adim + rind)),
                                    1 + 1/(2*(AA*R1_adim + rind)),
                                    np.ones(NG1+k, dtype='complex')])
        # Outer
        rows = np.append(rows, [np.arange(NG2+k, dtype='int')+(N1+1)*(M+1)+NG1+ringSubs+N2*k,
                                np.arange(NG2+k, dtype='int')+(N1+1)*(M+1)+NG1+ringSubs+N2*k,
                                np.arange(NG2+k, dtype='int')+(N1+1)*(M+1)+NG1+ringSubs+N2*k,
                                np.arange(NG2+k, dtype='int')+(N1+1)*(M+1)+NG1+ringSubs+N2*k,
                                np.arange(NG2+k, dtype='int')+(N1+1)*(M+1)+NG1+ringSubs+N2*k])
        cols = np.append(cols, [np.arange(NG2+k, dtype='int')+(N1+1)*M+stepSubs+NG1+ringSubs+N2*k+(k>0)*stepSubs,
                                np.arange(NG2+k, dtype='int')+(N1+1)*(M+1)+NG1+ringSubs+N2*k-1,
                                np.arange(NG2+k, dtype='int')+(N1+1)*(M+1)+NG1+ringSubs+N2*k,
                                np.arange(NG2+k, dtype='int')+(N1+1)*(M+1)+NG1+ringSubs+N2*k+1,
                                np.arange(NG2+k, dtype='int')+(N1+1)*(M+1)+NG1+ringSubs+N2*(k+1)+1])
        rind = np.arange(NG2+k, dtype='complex')+ringSubs+NG1-k
        coefs = np.append(coefs, [np.ones(NG2+k),
                                    1 - 1/(2*(AA*R1_adim + rind)),
                                    -1j*(Re1/(AA*AA)) - 4 - 1/((AA*R1_adim + rind)*(AA*R1_adim + rind)),
                                    1 + 1/(2*(AA*R1_adim + rind)),
                                    np.ones(NG2+k)])
    P2, Q2 = np.mgrid[ringSubs_2:M-1, 1:N2]
    rows = np.append(rows, [(M+1)*(N1+1)+(P2*(N2+1)+Q2).reshape(1, -1),
                            (M+1)*(N1+1)+(P2*(N2+1)+Q2).reshape(1, -1),
                            (M+1)*(N1+1)+(P2*(N2+1)+Q2).reshape(1, -1),
                            (M+1)*(N1+1)+(P2*(N2+1)+Q2).reshape(1, -1),
                            (M+1)*(N1+1)+(P2*(N2+1)+Q2).reshape(1, -1)])
    cols = np.append(cols, [(M+1)*(N1+1)+((P2-1)*(N2+1)+Q2).reshape(1, -1),
                            (M+1)*(N1+1)+(P2*(N2+1)+Q2).reshape(1, -1)-1,
                            (M+1)*(N1+1)+(P2*(N2+1)+Q2).reshape(1, -1),
                            (M+1)*(N1+1)+(P2*(N2+1)+Q2).reshape(1, -1)+1,
                            (M+1)*(N1+1)+((P2+1)*(N2+1)+Q2).reshape(1, -1)])
    rind = np.arange(N2-1, dtype='complex') + 1
    coefs = np.append(coefs, [np.ones((N2-1)*(M-ringSubs_2-1), dtype='complex'),
                                np.tile(1 - 1/(2*(AA*R1_adim + rind)), M-ringSubs_2-1),
                                np.tile(-1j*(Re1/(AA*AA)) - 4 - 1/((AA*R1_adim + rind)*(AA*R1_adim + rind)), M-ringSubs_2-1),
                                np.tile(1 + 1./(2*(AA*R1_adim + rind)), M-ringSubs_2-1),
                                np.ones((N2-1)*(M-ringSubs_2-1), dtype='complex')])
    ## Ring (no-slip)
    gDWR = np.array([], dtype='float64')
    indDWR = np.array([], dtype='int')
    for k in range(ringSubs_2+1):
        rows = np.append(rows, np.arange(1, 2*(k+1), dtype='int')+(M-ringSubs_2+k)*(N1+1)+NG1+stepSubs+ringSubs_2-k-1)
        cols = np.append(cols, np.arange(1, 2*(k+1), dtype='int')+(M-ringSubs_2+k)*(N1+1)+NG1+stepSubs+ringSubs_2-k-1)
        coefs = np.append(coefs, np.ones(2*(k+1)-1, dtype='float64'))
        gDWR = np.append(gDWR, R2_adim + (np.arange(2*(k+1)-1, dtype='float64')+NG1+stepSubs+ringSubs_2-k)/AA)
        indDWR = np.append(indDWR, np.arange(1, 2*(k+1), dtype='int')+(M-ringSubs_2+k)*(N1+1)+NG1+stepSubs+ringSubs_2-k-1)
        
    for k in range(ringSubs_2):
        rows = np.append(rows, np.arange(1, 2*(ringSubs_2-(k+1))+2, dtype='int')+NG1+1+(N1+1)*(M+1)+(N2+2)*k-1)
        cols = np.append(cols, np.arange(1, 2*(ringSubs_2-(k+1))+2, dtype='int')+NG1+1+(N1+1)*(M+1)+(N2+2)*k-1)
        coefs = np.append(coefs, np.ones((2*(ringSubs_2-(k+1))+1), dtype='float64'))
        gDWR = np.append(gDWR, R1_adim + (np.arange((2*(ringSubs_2-(k+1))+1), dtype='float64')+NG1+1+k)/AA)
        indDWR = np.append(indDWR, np.arange((2*(ringSubs_2-(k+1))+1), dtype='int')+NG1+1+(N1+1)*(M+1)+(N2+2)*k)

    ## FILL A
    A = sparse.csc_matrix((coefs, (rows, cols)), shape=(dim, dim))
    ## FILL b with no-slip boundary conditions over the DWR surfaces
    b = sparse.csc_matrix((gDWR, (indDWR, np.zeros(len(indDWR), dtype='int'))), shape=(dim, 1))
    # Solving the linear system of equations
    g = sparse.linalg.spsolve(A, b)

    # Torque calculation
    # indexes
    gpM = np.zeros((ringSubs + 5, ringSubs + 5), dtype='int')
    for k in range(ringSubs_2+3):
        gpM[k, :] = np.arange(5+ringSubs, dtype='int')+NG1+stepSubs+(M-ringSubs_2-2+k)*(N1+1)-2
    for k in range(ringSubs_2+2):
        gpM[k+ringSubs_2+3, :] = np.arange(5+ringSubs, dtype='int')+NG1+(N1+1)*(M+1)+(N2+1)*k-2
    gM = g[gpM]

    # Quadrants
    gp1 = np.zeros((3, ringSubs_2+1), dtype='complex')
    gp2 = np.zeros((3, ringSubs_2+1), dtype='complex')
    gp3 = np.zeros((3, ringSubs_2+1), dtype='complex')
    gp4 = np.zeros((3, ringSubs_2+1), dtype='complex')
    
    ind = np.arange(np.floor(((ringSubs+2)/2+1)/2), dtype='int')
    if np.remainder((ringSubs+2)/2+1, 2) != 0:
        ind = np.append(ind, np.floor(((ringSubs+2)/2+1)/2)-1)
    ind = np.append(ind, np.arange(np.floor(((ringSubs+2)/2+1)/2)-1, 0, -1, dtype='int')-1)

    for k in range(ringSubs_2 + 1):
        d1 = np.diag(gM, -ringSubs_2 + k*2)
        gp1[:, k] = d1[ind[k]:ind[k]+3]
        d2 = np.diag(np.rot90(gM), -ringSubs_2 + k*2)
        gp2[:, k] = d2[ind[k]:ind[k]+3]
        d3 = np.diag(np.rot90(gM, 2), -ringSubs_2 + k*2)
        gp3[:, k] = d3[ind[k]:ind[k]+3]
        d4 = np.diag(np.rot90(gM, 3), -ringSubs_2 + k*2)
        gp4[:, k] = d4[ind[k]:ind[k]+3]
        
    r1 = R2_adim + np.arange(stepSubs+NG1, stepSubs+0.5*ringSubs+NG1+1)*Sres_adim
    T1_fo = (1/(np.sqrt(2)))*np.trapz(r1*r1*(gp1[1, :] - gp1[2, :])) # first order
    T1_so = (1/(2*np.sqrt(2)))*np.trapz(r1*r1*(-3*gp1[2, :] + 4*gp1[1, :] - gp1[0, :])) # second order
    T1 = np.array([T1_fo, T1_so])
    
    r2 = R2_adim + np.arange(stepSubs+0.5*ringSubs+NG1, stepSubs+ringSubs+NG1+1)*Sres_adim
    T2_fo = (1/(np.sqrt(2)))*np.trapz(r2*r2*(gp2[1, :] - gp2[2, :]))
    T2_so = (1/(2*np.sqrt(2)))*np.trapz(r2*r2*(-3*gp2[2, :] + 4*gp2[1, :] - gp2[0, :]))
    T2 = np.array([T2_fo, T2_so])

    r3 = np.flip(r2)
    T3_fo = (1/(np.sqrt(2)))*np.trapz(r3*r3*(gp3[1, :] - gp3[2, :]))
    T3_so = (1/(2*np.sqrt(2)))*np.trapz(r3*r3*(-3*gp3[2, :] + 4*gp3[1, :] - gp3[0, :]))
    T3 = np.array([T3_fo, T3_so])
    
    r4 = np.flip(r1)
    T4_fo = (1/(np.sqrt(2)))*np.trapz(r4*r4*(gp4[1, :] - gp4[2, :]))
    T4_so = (1/(2*np.sqrt(2)))*np.trapz(r4*r4*(-3*gp4[2, :] + 4*gp4[1, :] - gp4[0, :]))
    T4 = np.array([T4_fo, T4_so])
    
    gsin = gM[ringSubs_2+2, 0:3]
    rin = R2_adim + np.arange(stepSubs+NG1-2, stepSubs+NG1+1)*Sres_adim
    g_r_in = gsin/rin
    Ts_in_fo = R5_adim*R5_adim*R5_adim*(1/(Sres_adim))*(g_r_in[1] - g_r_in[2])
    Ts_in_so = R5_adim*R5_adim*R5_adim*(1/(2*Sres_adim))*(-3*g_r_in[2] + 4*g_r_in[1] - g_r_in[0]);
    Ts_in = np.array([Ts_in_fo, Ts_in_so])
    
    gsout = gM[ringSubs_2+2, ringSubs+2:ringSubs+5]
    rout = R2_adim + np.arange(stepSubs+NG1+ringSubs, stepSubs+NG1+ringSubs+3)*Sres_adim
    g_r_out = gsout/rout
    Ts_out_fo = (1/(Sres_adim))*(g_r_out[1] - g_r_out[0])
    Ts_out_so = (1/(2*Sres_adim))*(-3*g_r_out[0] + 4*g_r_out[1] - g_r_out[2])
    Ts_out = np.array([Ts_out_fo, Ts_out_so])

    return g, T1, T2, T3, T4, Ts_in, Ts_out

def postprocessingDWR_ll(geom, mesh, Bulk, iteParams, IO):
# This code calculates the rheological properties from amplitude ratio
# and phase lag data of a rotational rheometer with a DWR fixture

# INPUTS (all units are in the International System)
# geom: geometry parameters
# mesh: mesh parameters
# bulk: bulk phases physical parameters
# iteParams: iterative scheme parameters
# iteMax: maximum number of iterations
# IO: Input/Output data parameters

# OUTPUTS(optional)
# resultStruct: strut with results
# timeElapsedTotal: time elapsed in processing experiental files

    start_timeTot = time.time()
    resultStruct = struct()

    expFilenames = []
    if sys.platform == 'win32':
        for name  in glob.glob(IO.inputFilepath + '\\*_exp.txt'):
            expFilenames.append(name)
    else:
        for name  in glob.glob(inputFilepath + '/*_exp.txt'):
            expFilenames.append(name)
            
    # Initializing optional output data
    iterationsTimesData = []
    iterationsData = []
    bouData = []
    FreqData = []
    AmpData = []
    TorqData = []
    etasData = []
    etasData_linear = []
    GData = []
    GData_linear = []
    ARcalcData = []
    deltaARcalcData = []
    # Displaying the number of input data files to process
    print('Experiment files: {} \n'.format(len(expFilenames)))

    # for-end looping on each input data file (*_exp.txt file)
    for ite in range(len(expFilenames)):
        print('Analyzing file {}...\n'.format(ite+1))
        if sys.platform == 'win32':
            expFileData = np.genfromtxt(expFilenames[ite].split('\\')[-1], delimiter='\t', dtype=None)
        else:
            expFileData = np.genfromtxt(expFilenames[ite].split('/')[-1], delimiter='\t', dtype=None)
            
        if len(expFileData.shape) == 1:
            expFileData = np.array([expFileData])    
        # Prompting the number of lines to process in the file
        expFileData = expFileData.T
        filasexpFiledata = int(expFileData[0].shape[0])        
        print('Lines: {} \n'.format(filasexpFiledata))

        # for-end looping on each line under analysis
        timeElapsedIT = np.zeros(filasexpFiledata, 'float64')
        Bou_final = np.zeros(filasexpFiledata, 'complex')
        lambda_final = np.zeros(filasexpFiledata, 'int')
        ARcalc_final = np.zeros(filasexpFiledata, 'complex')
        errorAR_final = np.zeros(filasexpFiledata, 'float64')
        delta_AR_final = np.zeros(filasexpFiledata, 'float64')
        freq = np.zeros(filasexpFiledata, 'float64')
        amp = np.zeros(filasexpFiledata, 'float64')
        Torq = np.zeros(filasexpFiledata, 'float64')
        G_linear = np.zeros((filasexpFiledata, 2), 'complex')
        etas_linear = np.zeros((filasexpFiledata, 2), 'complex')        

        for lin in range(filasexpFiledata):
            freq[lin] = expFileData[IO.colIndexFreq][lin]
            amp[lin] = expFileData[IO.colIndexAmp][lin]
            Torq[lin] = expFileData[IO.colIndexTorq][lin]
            omegarad = 2*np.pi*freq[lin]
            # Displaying the number of the line under analysis
            print('Analyzing line {} from file {}...\n'.format(lin+1, ite+1))
            Re1 = (Bulk.rho_bulk1*omegarad*geom.R6*geom.R6)/Bulk.eta_bulk1# Reynolds number
            Re2 = (Bulk.rho_bulk2*omegarad*geom.R6*geom.R6)/Bulk.eta_bulk2# Reynolds number
            Y = Bulk.eta_bulk1/Bulk.eta_bulk2
            ARexp = expFileData[IO.colIndexAR][lin]*(np.cos(expFileData[IO.colIndexDelta][lin]) + 1j*np.sin(expFileData[IO.colIndexDelta][lin]))
            # Intializing variables...
            lmbd = 1
            Bou = []
            NN = []
            ARcalc = []
            errorAR = []            
            Bou.append(0.0)
            start_timeIt = time.time()
            # Solving the Navier-Stokes equation
            NN.append(Bou[-1]*(Bulk.eta_bulk1 + Bulk.eta_bulk2)/Bulk.eta_bulk1)
            R5_adim = geom.R5/geom.R6
            g, T1, T2, T3, T4, Ts_in, Ts_out = solve_NS_DWR_ll(Re1, Re2, NN[-1], Y, geom.H, geom.G1, geom.G2, geom.R1, geom.R2, geom.R6, geom.ringW, geom.stepW, mesh.ringSubs, mesh.upperBC, R5_adim)
            # Calculating AR
            if mesh.DOrder == 1:
                T1 = T1[0]; T2 = T2[0]; T3 = T3[0]; T4 = T4[0]; Ts_in = Ts_in[0]; Ts_out = Ts_out[0]
            else:
                T1 = T1[1]; T2 = T2[1]; T3 = T3[1]; T4 = T4[1]; Ts_in = Ts_in[1]; Ts_out = Ts_out[1]
            Cll = 1j*omegarad*2*np.pi*geom.R6**3
            Tdown = Cll*Bulk.eta_bulk1*(T1 + T2)
            Tup = Cll*Bulk.eta_bulk2*(T3 + T4)
            Ts = Cll*NN[-1]*Bulk.eta_bulk1*(Ts_in + Ts_out)
            if geom.ICorrected:
                ARcalc.append(-Tdown - Tup - Ts)
            else:
                ARcalc.append(-Tdown - Tup - Ts - geom.inertia*omegarad**2)

            # Calculating first value of the relative error 
            errorAR.append(np.abs((ARcalc[-1])/(ARexp)-1))
                        
            # while loop performing the successive iterations     
            while lmbd < iteParams.iteMax and errorAR[-1] > iteParams.tolMin: 
                # print(lmbd)
                lmbd += 1
                if geom.ICorrected:
                    Bou.append((-ARexp - Tdown - Tup)/(Cll*(Bulk.eta_bulk1 + Bulk.eta_bulk2)*(Ts_in + Ts_out)))
                else:
                    Bou.append((-ARexp - Tdown - Tup - geom.inertia*omegarad**2)/(Cll*(Bulk.eta_bulk1 + Bulk.eta_bulk2)*(Ts_in + Ts_out)))          
                # Solving the Navier-Stokes equation
                NN.append(Bou[-1]*(Bulk.eta_bulk1 + Bulk.eta_bulk2)/Bulk.eta_bulk1)
                g, T1, T2, T3, T4, Ts_in, Ts_out = solve_NS_DWR_ll(Re1, Re2, NN[-1], Y, geom.H, geom.G1, geom.G2, geom.R1, geom.R2, geom.R6, geom.ringW, geom.stepW, mesh.ringSubs, mesh.upperBC, R5_adim)
                # Calculating AR
                if mesh.DOrder == 1:
                    T1 = T1[0]; T2 = T2[0]; T3 = T3[0]; T4 = T4[0]; Ts_in = Ts_in[0]; Ts_out = Ts_out[0]
                else:
                    T1 = T1[1]; T2 = T2[1]; T3 = T3[1]; T4 = T4[1]; Ts_in = Ts_in[1]; Ts_out = Ts_out[1]
                Tdown = Cll*Bulk.eta_bulk1*(T1 + T2)
                Tup = Cll*Bulk.eta_bulk2*(T3 + T4)
                Ts = Cll*NN[-1]*Bulk.eta_bulk1*(Ts_in + Ts_out)
                if geom.ICorrected:
                    ARcalc.append(-Tdown - Tup - Ts)
                else:
                    ARcalc.append(-Tdown - Tup - Ts - geom.inertia*omegarad**2)
                # Calculating first value of the relative error 
                errorAR.append(np.abs((ARcalc[-1])/(ARexp)-1))
                
            # Exit of the iterative process
            if lmbd == iteParams.iteMax:
                print('Limit reached \n')# Maximum # of iterations reached
            else:
                print('OK convergence at {} iterations\n'.format(lmbd))# Iterations have converged!!!
                
            timeElapsedIT[lin] = time.time() - start_timeIt
            # Displaying the time used in the iterative process for each line
            print('Iterative process time = {} s\n'.format(timeElapsedIT[lin]))
            Bou_final[lin] = Bou[-1]
            ARcalc_final[lin] = ARcalc[-1]
            delta_AR_final[lin] = np.arctan(np.imag(ARcalc[-1])/np.real(ARcalc[-1]))
            if delta_AR_final[lin] < 0:
                delta_AR_final[lin] = delta_AR_final[lin] + np.pi
            errorAR_final[lin] = errorAR[-1]
            lambda_final[lin] = lmbd
            G_linear[lin][0] = ARexp/(4*np.pi*((geom.R5**2*geom.R1**2/(geom.R5**2 - geom.R1**2)) + (geom.R6**2*geom.R3**2/(geom.R3**2 - geom.R6**2))))
            etas_linear[lin][0] = G_linear[lin, 0]/(1j*omegarad)
            G_linear[lin][1] = (ARexp - ARcalc[0])/(4*np.pi*((geom.R5**2*geom.R1**2/(geom.R5**2 - geom.R1**2)) + (geom.R6**2*geom.R3**2/(geom.R3**2 - geom.R6**2))))
            etas_linear[lin][1] = G_linear[lin, 1]/(1j*omegarad)

        # Calculating variables depending on the N number
        eta_s_final = Bou_final*geom.R6*(Bulk.eta_bulk1 + Bulk.eta_bulk2)# converged viscoelasticity
        G_complex = 1j*Bou_final*2*np.pi*freq*geom.R6*(Bulk.eta_bulk1 + Bulk.eta_bulk2)# converged dynamic surface moduli    
        
        # Exporting results to the output data file
        # aa = np.array(expFileData)
        # bb = [tup[0] for tup in aa]
        # cc = [tup[IO.colIndexFreq] for tup in bb]
        
        results = (np.array([expFileData[IO.colIndexFreq], np.real(G_complex), np.imag(G_complex), np.real(eta_s_final),
                                np.imag(eta_s_final), np.real(Bou_final), np.imag(Bou_final), np.absolute(ARcalc_final),
                                delta_AR_final, timeElapsedIT, lambda_final]))
        results = results.T
        if sys.platform == 'win32':
            np.savetxt(IO.outputFilepath+'\\'+expFilenames[ite].split('\\')[-1].replace('exp','out'), results, fmt='%.14f', delimiter='\t')
        else:
            np.savetxt(IO.outputFilepath+'/'+expFilenames[ite].split('/')[-1].replace('exp','out'), results, fmt='%.14f', delimiter='\t')
        
        # Optional structure output data
        FreqData.append(freq)
        AmpData.append(amp)
        TorqData.append(Torq)
        GData.append(G_complex)
        GData_linear.append(G_linear)
        etasData.append(eta_s_final)
        etasData_linear.append(etas_linear)
        bouData.append(Bou_final)
        ARcalcData.append(np.absolute(ARcalc_final))
        deltaARcalcData.append(delta_AR_final)
        iterationsData.append(lambda_final)
    resultStruct.FreqData = FreqData
    resultStruct.AmpData = AmpData
    resultStruct.TorqData = TorqData
    resultStruct.GData = GData
    resultStruct.GData_linear = GData_linear
    resultStruct.etasData = etasData
    resultStruct.etasData_linear = etasData_linear
    resultStruct.bouData = bouData
    resultStruct.ARcalcData = ARcalcData
    resultStruct.deltaARcalcData = deltaARcalcData
    resultStruct.iterationsTimesData = iterationsTimesData
    resultStruct.iterationsData = iterationsData
    resultStruct.geom = geom
    resultStruct.mesh = mesh
    resultStruct.Bulk = Bulk
    resultStruct.iteParams = iteParams
    timeElapsedTotal = time.time() - start_timeTot
    print('Total postprocessing program time = {}\n'.format(timeElapsedTotal))# Iterations have converged!!!

    return resultStruct