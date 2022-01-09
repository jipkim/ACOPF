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

#%% pyomo] start
M = pyo.ConcreteModel() 
M.theta = pyo.Var( busset, bounds=( - np.pi, + np.pi ) ) 
M.pg = pyo.Var( genset, bounds = pg_bound )
M.p_ft = pyo.Var( lineset ) 
M.p_tf = pyo.Var( lineset )
#%% constrains
for b in busset:
    if buses[b].btype == 3:
        M.theta[b].fix(0)

def p_ft_def( M, l ):
    return M.p_ft[l] == float(yft_i[l]) * (M.theta[lines[l].fbus] - M.theta[lines[l].tbus])
def p_tf_def(M, l):
    return M.p_tf[l] == float(ytf_i[l]) * (M.theta[lines[l].tbus] - M.theta[lines[l].fbus])

M.Con_p_ft_def = pyo.Constraint( lineset, rule=p_ft_def )
M.Con_p_tf_def = pyo.Constraint( lineset, rule=p_tf_def )

def Pbalance( M, b ):
    return sum(M.p_ft[l] for l in buses[b].outline) \
        + sum(M.p_tf[l] for l in buses[b].inline) \
        - sum(M.pg[g] for g in B_gn[b]) \
        + buses[b].Pd \
        + buses[b].Gs * 1.0 ** 2 \
        == 0
M.Con_Pbalance = pyo.Constraint( busset, rule=Pbalance )



def p_ft_ub( M, l ):
    return M.p_ft[l] <= + lines[l].u
def p_ft_lb( M, l ):
    return M.p_ft[l] >= - lines[l].u
def p_tf_ub(M, l):
    return M.p_tf[l] <= + lines[l].u
def p_tf_lb(M, l):
    return M.p_tf[l] >= - lines[l].u

lineset_thermal = [];
for l in lineset:
    if lines[l].u != 0:
        lineset_thermal.append(l)

M.Con_p_ft_ub = pyo.Constraint(lineset_thermal, rule=p_ft_ub)
M.Con_p_ft_lb = pyo.Constraint(lineset_thermal, rule=p_ft_lb)
M.Con_p_tf_ub = pyo.Constraint(lineset_thermal, rule=p_tf_ub)
M.Con_p_tf_lb = pyo.Constraint(lineset_thermal, rule=p_tf_lb)


def ang_diff_ub(M, l):
    # return M.theta[lines[l].fbus] - M.theta[lines[l].tbus] <= lines[l].angmax
    return M.theta[lines[l].fbus] - M.theta[lines[l].tbus] <= min(lines[l].angmax, radians(60))
def ang_diff_lb(M, l):
    # return M.theta[lines[l].fbus] - M.theta[lines[l].tbus] >= lines[l].angmin
    return M.theta[lines[l].fbus] - M.theta[lines[l].tbus] >= max(lines[l].angmin, -radians(60))

M.Con_ang_diff_ub = pyo.Constraint(lineset, rule=ang_diff_ub)
M.Con_ang_diff_lb = pyo.Constraint(lineset, rule=ang_diff_lb)


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
