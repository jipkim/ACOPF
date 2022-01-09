#%%
import numpy as np
import cvxpy as cp

# from dataclasses import dataclass, field
# from typing import List
from math import radians

#%%
# from src.struct_network import Bus, Line, Generator
from src.readnetwork import readnetwork
from src.ybus import ybus
#%%


testsystem = "case_ieee30"

#%%
buses, lines, generators, datamat = readnetwork(testsystem)

nline = len(lines)
nbus = len(buses)
ngen = len(generators)
lineset = np.arange( 0, nline )
busset = np.arange( 0, nbus )
genset = np.arange(0, ngen )
baseMVA = datamat["baseMVA"][0, 0][0, 0]

B_g = []
for g in genset:
    B_g.append( generators[g].location )

B_gn = []
for i in busset:
    B_gn.append( [] )

for g in genset:
    B_gn[ generators[g].location ].append( int(g) )


#%%
Ybus, yff, yft, ytf, ytt = ybus( buses, lines )

yff_r = np.real(yff)
yff_i = np.imag(yff)
ytt_r = np.real(ytt)
ytt_i = np.imag(ytt)
yft_r = np.real(yft)
yft_i = np.imag(yft)
ytf_r = np.real(ytf)
ytf_i = np.imag(ytf)



#%% cvx start

theta = cp.Variable( len(busset))
pg = cp.Variable( len(genset) )
p_ft = cp.Variable( len(lineset) )
p_tf = cp.Variable( len(lineset) )

constraints = []
for b in busset:
    constraints += [ - np.pi <= theta[b] ]
    constraints += [ theta[b] <= + np.pi ] 
    if buses[b].btype == 3:
        constraints += [ theta[b] == 0]
    
for g in genset:
    constraints += [ generators[g].Pmin <= pg[g] ]
    constraints += [ pg[g] <= generators[g].Pmax ]

for l in lineset:
    constraints += [ 
        p_ft[l] == float(yft_i[l]) * (theta[lines[l].fbus]-theta[lines[l].tbus]) 
            ]
    constraints += [ 
        p_tf[l] == float(ytf_i[l]) * (theta[lines[l].tbus]-theta[lines[l].fbus])
            ]


for b in busset:
    constraints += [ 
        sum(p_ft[l] for l in buses[b].outline) 
        + sum(p_tf[l] for l in buses[b].inline) 
        - sum(pg[g] for g in B_gn[b]) 
        + buses[b].Pd 
        + buses[b].Gs * 1.0 ** 2 
        == 0
        ]

for l in lineset:
    if lines[l].u != 0:
        constraints += [ p_ft[l] <= + lines[l].u ]
        constraints += [ p_ft[l] >= - lines[l].u ]
        constraints += [ p_tf[l] <= + lines[l].u ]
        constraints += [ p_tf[l] >= - lines[l].u ]    
    constraints += [ theta[lines[l].fbus] - theta[lines[l].tbus] <= lines[l].angmax ]
    constraints += [ theta[lines[l].fbus] - theta[lines[l].tbus] >= lines[l].angmin ]
    # constraints += [ theta[lines[l].fbus] - theta[lines[l].tbus] <= min(lines[l].angmax, radians(60)) ]
    # constraints += [ theta[lines[l].fbus] - theta[lines[l].tbus] >= max(lines[l].angmin, -radians(60)) ]

prob = cp.Problem(cp.Minimize(
    sum(
        (generators[g].cost[0] * (pg[g]*baseMVA)**2
        + generators[g].cost[1] * (pg[g]*baseMVA) 
        + generators[g].cost[2]) 
            for g in genset) ), constraints)


#%% 
# print(cp.installed_solvers())
prob.solve( solver=cp.GUROBI , verbose = True )
# prob.solve( solver=cp.MOSEK , verbose = True )
# prob.solve( verbose = True )
print("Termination status =", prob.status)
print("The optimal value is", prob.value)
print("test system =", testsystem)
print("solve time =", prob.solver_stats.solve_time)
