using JuMP, Ipopt
# using Gurobi
using SparseArrays # ybus
# using DataFrames # for comparison 


include(string(pwd(), "/src/struct_network.jl"))
include(string(pwd(), "/src/readnetwork.jl"))
include(string(pwd(), "/src/ybus.jl"))
#----------------------------------------------------------------------------------#
solver = "Ipopt"
# solver = "Gurobi"

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


############
Ybus, yff, yft, ytf, ytt = ybus( buses, lines )

yff_r = real.(yff)
yff_i = imag.(yff)
ytt_r = real.(ytt)
ytt_i = imag.(ytt)
yft_r = real.(yft)
yft_i = imag.(yft)
ytf_r = real.(ytf)
ytf_i = imag.(ytf)



#---------------------------------------------------------------------------------#
if solver == "Gurobi"
   m = Model(Gurobi.Optimizer)
   # set_optimizer_attribute(m, "MIPGap", 1e-4)
   # set_optimizer_attribute(m, "Method", 5)
   # set_optimizer_attribute(m, "Presolve", 0) # max 2
   # set_optimizer_attribute(m, "NonConvex", 2)
#    set_optimizer_attribute(m, "NumericFocus", 3) # max 3
   # set_optimizer_attribute(m, "Threads", 6) 
elseif solver == "Mosek"
   m = Model(Mosek.Optimizer)
elseif solver == "KNITRO"
   m = Model(KNITRO.Optimizer)   
elseif solver == "Ipopt"
   m = Model(Ipopt.Optimizer)   
end

@variables(m, begin
   - π <= θ[b in busset] <= π
   generators[g].Pmin <= pg[g in genset] <= generators[g].Pmax  # real power output of generator
end)

for b in busset 
    if buses[b].btype == 3
        fix(θ[b], 0; force=true) # θ[1] = 0 at root node
    end
end



@expression(m, p_ft[l in lineset], + yft_i[l]*(θ[lines[l].fbus]-θ[lines[l].tbus]))
@expression(m, p_tf[l in lineset], + ytf_i[l]*(θ[lines[l].tbus]-θ[lines[l].fbus]))
    
@constraint(m, PBalance[b in busset],     
    + sum(p_ft[l] for l in buses[b].outline)  
    + sum(p_tf[l] for l in buses[b].inline)     
    - sum(pg[g] for g in B_gn[b]) 
    + buses[b].Pd 
    + buses[b].Gs*1.0^2
        == 0)

@constraint(m, LineCap_ub[l in lineset; lines[l].u != 0], p_ft[l] <= + lines[l].u)
@constraint(m, LineCap_lb[l in lineset; lines[l].u != 0], p_ft[l] >= - lines[l].u)
@constraint(m, LineCapBW_ub[l in lineset; lines[l].u != 0], p_tf[l] <= + lines[l].u)
@constraint(m, LineCapBW_lb[l in lineset; lines[l].u != 0], p_tf[l] >= - lines[l].u)

@constraint(m, AngDiffmax[l in lineset], θ[lines[l].fbus] - θ[lines[l].tbus] <= min(lines[l].angmax, deg2rad(60)))
@constraint(m, AngDiffmin[l in lineset], θ[lines[l].fbus] - θ[lines[l].tbus] >= max(lines[l].angmin, -deg2rad(60)))


@expression(m, OCgen[g in genset], 
    generators[g].cost[1]*(pg[g]*datamat["baseMVA"])^2 
    + generators[g].cost[2]*(pg[g]*datamat["baseMVA"]) 
    + generators[g].cost[3])
@expression(m, OC, sum(OCgen[g] for g in genset))

@objective(m, Min, OC)


optimize!(m)
@show termination_status(m)
@show objective_value(m)


