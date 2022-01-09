#%%
import pyomo.environ as pyo
import numpy as np
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

#%%
def pg_bound(m, g):
    return(generators[g].Pmin, generators[g].Pmax)
def qg_bound(m, g):
    return(generators[g].Qmin, generators[g].Qmax)
def u_bound(m, b):
    return(buses[b].Vmin**2, buses[b].Vmax**2)
#%% pyomo] start
M = pyo.ConcreteModel() 
M.u = pyo.Var( busset, bounds = u_bound )
M.a = pyo.Var( lineset, domain = pyo.NonNegativeReals )
M.pg = pyo.Var( genset, bounds = pg_bound )
M.qg = pyo.Var( genset, bounds = qg_bound )
M.fp = pyo.Var(lineset)
M.fq = pyo.Var(lineset)


#%% constrains
# for b in busset:
#     if buses[b].btype == 3:        
#         M.v[b].fix(1)

def BetweenNodes( M, l ):
    return M.u[lines[l].tbus] \
        - 2*(lines[l].r*M.fp[l] + lines[l].x*M.fq[l]) \
        + M.a[l]*(lines[l].r**2 + lines[l].x**2) \
        == M.u[lines[l].fbus]

M.Con_BetweenNodes = pyo.Constraint( lineset, rule=BetweenNodes )

def PBalance( M, b ):
    return sum(M.fp[l] for l in buses[b].inline) \
        - sum(M.fp[l] - lines[l].r*M.a[l] for l in buses[b].outline) \
        - sum(M.pg[g] for g in B_gn[b]) \
        + buses[b].Pd \
        + buses[b].Gs*M.u[b] \
        == 0

def QBalance( M, b ):
    return sum(M.fq[l] for l in buses[b].inline) \
        - sum(M.fq[l] - lines[l].x*M.a[l] for l in buses[b].outline) \
        - sum(M.qg[g] for g in B_gn[b]) \
        + buses[b].Qd \
        - buses[b].Bs*M.u[b] \
        == 0

M.Con_PBalance = pyo.Constraint( busset, rule=PBalance )
M.Con_QBalance = pyo.Constraint( busset, rule=QBalance )



def LineCap_FW(M, l):
    return M.fp[l]**2 + M.fq[l]**2 <= + lines[l].u**2
def LineCap_BW(M, l):
    return (M.fp[l] - M.a[l]*lines[l].r)**2 + (M.fq[l] - M.a[l]*lines[l].x)**2 <= (lines[l].u)**2

lineset_thermal = [];
for l in lineset:
    if lines[l].u != 0:
        lineset_thermal.append(l)
M.Con_LineCap_FW = pyo.Constraint(lineset_thermal, rule=LineCap_FW)
M.Con_LineCap_BW = pyo.Constraint(lineset_thermal, rule=LineCap_BW)

def SOCP(M, l):
    return (2*M.fp[l])**2 + (2*M.fq[l])**2 + (M.u[lines[l].tbus] - M.a[l])**2 <= (M.u[lines[l].tbus] + M.a[l])**2
M.Con_SOCP = pyo.Constraint( lineset, rule=SOCP )


def obj_func(M):
    return sum(
        (generators[g].cost[0] * (M.pg[g]*baseMVA)**2
         + generators[g].cost[1] * (M.pg[g]*baseMVA)
         + generators[g].cost[2])
        for g in genset)
M.Obj = pyo.Objective(rule=obj_func, sense=pyo.minimize)

#%%

solvername = 'gurobi'
# solvername = 'mosek'
# solvername = 'ipopt'
# solvername = 'knitroampl'
solver = pyo.SolverFactory(solvername)
# solver.options["threads"] = 4
# solver.options["tee"] = True
# results = solver.solve(M)
if solvername == 'ipopt':
    solver.options["max_iter"] =  1_000_000
    solver.options["linear_solver"] =  "ma86"
    solver.options["print_level"] =  5  # default = 5 (0~12)
elif solvername == 'gurobi':
    # solver.options["MIPGap"] = 1e-6
    # solver.options["Method"] = 5
    # solver.options["Presolve"] = 0 # max 2
    # solver.options["NonConvex"] = 2
    # solver.options["NumericFocus"] = 3 # max 3
    solver.options["Threads"] = 8
elif solvername == 'knitroampl':
    # solver.options["option_file"] = "knitro.opt"
    solver.options["ms_enable"] = 0
    # solver.options["ms_maxsolves"] = 300
    # solver.options["ms_terminate"] = 1 # 0 = maxsolve | 1 = first local optimum
    # solver.options["par_blasnumthreads"] = 10
    # solver.options["maxit"] = 100000
    # solver.options["algorithm"] = 1
    # solver.options["convex"] = 1

results = solver.solve( M, tee=True)

print("Solver =", solvername)
print("Termination status =", results.solver.termination_condition)
print("The optimal value is", pyo.value(M.Obj))
print("test system =", testsystem)
