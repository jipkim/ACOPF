#%%
import cvxpy as cp
import numpy as np
from math import radians

#%%
# from src.struct_network import Bus, Line, Generator
from src.readnetwork import readnetwork
from src.ybus import ybus
from src.sdp_ymat import sdp_ymat


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
Yk, Yk_, Mk, Ylineft, Ylinetf, Y_lineft, Y_linetf, YL, YL_ = sdp_ymat(
    lines, Ybus)


gcostdata = datamat["gencost"][0, 0]

if np.size(gcostdata, 2-1) == 7:
    c_k2 = gcostdata[:, 5-1] * baseMVA ** 2
    c_k1 = gcostdata[:, 6-1] * baseMVA
    c_k0 = gcostdata[:, 7-1].T
elif np.size(gcostdata, 2-1) == 6:
    c_k2 = np.zeros(np.size(gcostdata, 1-1), 1)
    c_k1 = gcostdata[:, 5-1] * baseMVA
    c_k0 = gcostdata[:, 6-1].T

#%%
# p5 right column: a small resistance (1e-5) is added to each transformer ...
# ∴ connected graph
# for l in lineset:
#     if lines[l].r == 0:
#         lines[l].r == 1e-5;
#%% cvx start

# W = cp.Variable((2*len(busset), 2*len(busset)), symmetric=True)
W = cp.Variable((2*len(busset), 2*len(busset)), PSD=True)
alpha_k = cp.Variable(len(genset))

constraints = []
for b in busset:
    if buses[b].btype == 3:
        constraints += [W[:,nbus+b] == 0]
        constraints += [W[nbus+b,:] == 0]
        slack_bus = b;
        print("slack bus:", slack_bus)

for b in busset:
    # Equations (4a)
    constraints += [sum(generators[g].Pmin for g in B_gn[b]) - buses[b].Pd <= cp.trace(Yk(b) @ W)]
    constraints += [sum(generators[g].Pmax for g in B_gn[b]) - buses[b].Pd >= cp.trace(Yk(b) @ W)]
    # Equations (4b)
    constraints += [sum(generators[g].Qmin for g in B_gn[b]) - buses[b].Qd <= cp.trace(Yk_(b) @ W)]
    constraints += [sum(generators[g].Qmax for g in B_gn[b]) - buses[b].Qd >= cp.trace(Yk_(b) @ W)]
    # Equations (4c)
    constraints += [(buses[b].Vmin)**2 <= cp.trace(Mk(b) @ W)]
    constraints += [(buses[b].Vmax)**2 >= cp.trace(Mk(b) @ W)]


for l in lineset:
    if lines[l].u != 0:
        # Equations (4e)
        constraints += [- lines[l].u <= cp.trace(Ylineft(l) @ W)]
        constraints += [+ lines[l].u >= cp.trace(Ylineft(l) @ W)]
        # Equations (5): 3x3 matrix
        constraints += [
            cp.bmat([
                [-(lines[l].u)**2, cp.trace(Ylineft(l) @ W), cp.trace(Y_lineft(l) @ W)],
                [cp.trace(Ylineft(l) @ W), -1, 0],
                [cp.trace(Y_lineft(l) @ W), 0, -1]
            ]) << 0]

for g in genset:
    # Equations (6): 2x2 matrix
    constraints += [
        cp.bmat([
            [
                c_k1[g] * cp.trace(Yk(B_g[g]) @ W) - alpha_k[g] + c_k0[g] + c_k1[g] * buses[B_g[g]].Pd,
                cp.sqrt(c_k2[g]) * cp.trace(Yk(B_g[g]) @ W) + cp.sqrt(c_k2[g]) * (buses[B_g[g]].Pd)
            ],
            [
                cp.sqrt(c_k2[g]) * cp.trace(Yk(B_g[g]) @ W) + cp.sqrt(c_k2[g]) * (buses[B_g[g]].Pd), 
                -1
            ]
        ]) << 0]
    

#%%
prob = cp.Problem(cp.Minimize(sum(alpha_k[g] for g in genset)), constraints)

#%%
# print(cp.installed_solvers())
# prob.solve(solver=cp.GUROBI, verbose=True)
prob.solve( solver=cp.MOSEK , verbose = True )
# prob.solve( verbose = True )
print("Termination status =", prob.status)
print("The optimal value is", prob.value)
print("test system =", testsystem)
print("solve time =", prob.solver_stats.solve_time)


#%%
v = [np.sqrt(np.trace(Mk(b) @ W.value)) for b in busset]
fp = [np.trace(Ylineft(l) @ W.value) for l in lineset]
fq = [np.trace(Y_lineft(l) @ W.value) for l in lineset]
pg = [np.trace(Yk(B_g[g]) * W.value) + buses[B_g[g]].Pd for g in genset]
qg = [np.trace(Yk_(B_g[g]) * W.value) + buses[B_g[g]].Qd for g in genset]


# for variable in prob.variables():
#     print("Variable %s: value %s" % (variable.name(), variable.value))
