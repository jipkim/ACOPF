using JuMP, Ipopt
using Gurobi
using DataFrames # for comparison 

include(string(pwd(), "/src/struct_network.jl"))
include(string(pwd(), "/src/readnetwork.jl"))
#----------------------------------------------------------------------------------#
solver = "Gurobi"
# solver = "Ipopt"


testsystem = "case_ieee30"

buses, lines, generators, datamat = readnetwork(testsystem)

lineset = 1:length(lines)
busset = 1:length(buses)
genset = 1:length(generators)

B_g = []
for g in genset
   push!(B_g, generators[g].location)
end

B_gn = Array{Array{Int64,1}}(undef, length(buses))
for ii in 1:length(buses)
   B_gn[ii] = Int64[]
end
for g in genset
   push!(B_gn[generators[g].location], g)
end

#---------------------------------------------------------------------------------#
if solver == "Gurobi"
   m = Model(Gurobi.Optimizer)
   #    set_optimizer_attribute(m, "MIPGap", 1e-6)
   # set_optimizer_attribute(m, "Method", 5)
   # set_optimizer_attribute(m, "Presolve", 0) # max 2
   # set_optimizer_attribute(m, "NonConvex", 2)
   #    set_optimizer_attribute(m, "NumericFocus", 3) # max 3
   # set_optimizer_attribute(m, "Threads", 6) 
   # set_optimizer_attribute(m, "QCPdual", 1)
elseif solver == "Ipopt"
   m = Model(Ipopt.Optimizer)
   set_optimizer_attribute(m, "max_iter", 1_000_000)
   # set_optimizer_attribute(m, "linear_solver", "ma86")
   # set_optimizer_attribute(m, "print_level", 5) # default = 5 (0~12)
end


@variables(m, begin
   u[b in busset], (start = 1.0) # squared voltage magnitude, |v|^2
   a[l in lineset] >= 0             # squared current magnitude, |I|^2
   fp[l in lineset] # real power flow
   fq[l in lineset] # reactive power flow
   pg[g in genset]  # real power output of generator
   qg[g in genset]  # reactive power output of generator
end)


@constraint(m, PBalance[b in busset],
   sum(fp[l] for l in buses[b].inline)
   -
   sum(fp[l] + lines[l].r * a[l] for l in buses[b].outline)
   +
   sum(pg[g] for g in B_gn[b])
   -
   buses[b].Pd
   -
   buses[b].Gs * u[b]
   ==
   0)

@constraint(m, QBalance[b in busset],
   sum(fq[l] for l in buses[b].inline)
   -
   sum(fq[l] + lines[l].x * a[l] for l in buses[b].outline)
   +
   sum(qg[g] for g in B_gn[b])
   -
   buses[b].Qd
   +
   buses[b].Bs * u[b]
   ==
   0)

@constraint(m, BetweenNodes_reversed[l in lineset],
   u[lines[l].tbus]
   + 2 * (lines[l].r * fp[l] + lines[l].x * fq[l])
   + a[l] * (lines[l].r^2 + lines[l].x^2)
   ==
   u[lines[l].fbus])

@constraint(m, LineCapFW[l in lineset; lines[l].u != 0], (-fp[l])^2 + (-fq[l])^2 <= (lines[l].u)^2) # if lines[l].u == 0, then unlimited line capacity | matpower format
@constraint(m, LineCapBW[l in lineset; lines[l].u != 0], (-fp[l] - a[l] * lines[l].r)^2 + (-fq[l] - a[l] * lines[l].x)^2 <= (lines[l].u)^2)
# @constraint(m, SOCP[l in lineset], ((fp[l])^2 + (fq[l])^2) <= u[lines[l].tbus] * a[l] )
@constraint(m, SOCP[l in lineset], (-2 * fp[l])^2 + (-2 * fq[l])^2 + (u[lines[l].tbus] - a[l])^2 <= (u[lines[l].tbus] + a[l])^2)


# @constraint(m, LineCapFW[l in lineset; lines[l].u != 0], [lines[l].u, -fp[l], -fq[l]] in SecondOrderCone() ) # if lines[l].u == 0, then unlimited line capacity | matpower format
# @constraint(m, LineCapBW[l in lineset; lines[l].u != 0], [lines[l].u, -fp[l] - a[l]*lines[l].r, -fq[l] - a[l]*lines[l].x] in SecondOrderCone())
# @constraint(m, SOCP[l in lineset], [u[lines[l].tbus] + a[l], -2*fp[l], -2*fq[l], u[lines[l].tbus] - a[l]] in SecondOrderCone())

@constraint(m, PGmin[g in genset], pg[g] >= generators[g].Pmin)
@constraint(m, PGmax[g in genset], pg[g] <= generators[g].Pmax)
@constraint(m, QGmin[g in genset], qg[g] >= generators[g].Qmin)
@constraint(m, QGmax[g in genset], qg[g] <= generators[g].Qmax)

@constraint(m, Vmax[b in busset], u[b] <= buses[b].Vmax^2) #upper limit constraint for voltage square
@constraint(m, Vmin[b in busset], u[b] >= buses[b].Vmin^2) #lower limit constraint for voltage square

@expression(m, OCgen[g in genset],
   generators[g].cost[1] * (pg[g] * datamat["baseMVA"])^2
   + generators[g].cost[2] * (pg[g] * datamat["baseMVA"])
   + generators[g].cost[3])
@expression(m, OC, sum(OCgen[g] for g in genset))

@objective(m, Min, OC)


optimize!(m)
@show termination_status(m)
@show objective_value(m)
println()

