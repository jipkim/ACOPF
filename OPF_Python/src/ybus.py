import numpy as np
from scipy import sparse

def ybus(buses, lines):
    nbus = len(buses)
    nline = len(lines)
    busset = np.arange(0, nbus)
    lineset = np.arange(0,nline)


    Ybus = np.zeros( [nbus, nbus], dtype=complex ) 
    yff = np.zeros( [ nline, 1], dtype=complex )
    yft = np.zeros( [ nline, 1], dtype=complex )
    ytf = np.zeros( [ nline, 1], dtype=complex )
    ytt = np.zeros( [ nline, 1], dtype=complex )
    tau = np.zeros( [ nline, 1], dtype=complex )
    for l in lineset:
        
        tau[l] = np.exp(+ 1j * lines[l].shft) if lines[l].tap == 0 else lines[l].tap * np.exp(+ 1j * lines[l].shft)

        ytt[l] = 1/(lines[l].r + 1j * lines[l].x) + 1j * lines[l].b / 2
        yff[l] = (1/(lines[l].r + 1j * lines[l].x) + 1j * lines[l].b / 2) / (tau[l] * np.conj(tau[l]))
        yft[l] = - (1/(lines[l].r + 1j * lines[l].x)) / np.conj(tau[l])
        ytf[l] = - (1/(lines[l].r + 1j * lines[l].x)) / tau[l]

        Ybus[lines[l].fbus, lines[l].tbus] = yft[l]
        Ybus[lines[l].tbus, lines[l].fbus] = ytf[l]
        Ybus[lines[l].fbus, lines[l].fbus] += yff[l]
        Ybus[lines[l].tbus, lines[l].tbus] += ytt[l]
    
    for b in busset:
        Ybus[b, b] += buses[b].Gs + 1j * buses[b].Bs
    

    return Ybus, yff, yft, ytf, ytt

