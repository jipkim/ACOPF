using JuMP, Ipopt
# using KNITRO
using SparseArrays # ybus
# using DataFrames # for comparison 



include(string(pwd(), "/src/struct_network.jl"))
include(string(pwd(), "/src/readnetwork.jl"))
include(string(pwd(), "/src/ybus.jl"))
#----------------------------------------------------------------------------------#
solver = "Ipopt"
# solver = "KNITRO"


testsystem = "case_ieee30"
buses, lines, generators, datamat = readnetwork(testsystem)
#----------------------------------------------------------------------------------#

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

Ybus, yff, yft, ytf, ytt = ybus(buses, lines)

yff_r = real.(yff)
yff_i = imag.(yff)
ytt_r = real.(ytt)
ytt_i = imag.(ytt)
yft_r = real.(yft)
yft_i = imag.(yft)
ytf_r = real.(ytf)
ytf_i = imag.(ytf)




#---------------------------------------------------------------------------------#
if solver == "Ipopt"
    m = Model(Ipopt.Optimizer)
    set_optimizer_attribute(m, "max_iter", 1_000_000)
    # set_optimizer_attribute(m, "linear_solver", "ma86")
    # set_optimizer_attribute(m, "print_level", 5) # default = 5 (0~12)
elseif solver == "KNITRO"
    m = Model(KNITRO.Optimizer)
    # set_optimizer_attribute(m, "option_file", "knitro.opt")
    #    set_optimizer_attribute(m, "ms_enable", 0)
    #    set_optimizer_attribute(m, "ms_maxsolves", 300)
    #    set_optimizer_attribute(m, "ms_terminate", 1) # 0 = maxsolve | 1 = first local optimum
    #    set_optimizer_attribute(m, "par_blasnumthreads", 10)
    #    set_optimizer_attribute(m, "maxit", 100000)
    # set_optimizer_attribute(m, "algorithm", 5)
    # set_optimizer_attribute(m, "convex", 1)
end


@variables(m, begin
    buses[b].Vmin <= v[b in busset] <= buses[b].Vmax, (start = 1.0)
    - π <= θ[b in busset] <= π
    generators[g].Pmin <= pg[g in genset] <= generators[g].Pmax  # real power output of generator
    generators[g].Qmin <= qg[g in genset] <= generators[g].Qmax # reactive power output of generator
end)


for b in busset
    if buses[b].btype == 3
        fix(θ[b], 0; force = true) # θ[1] = 0 at root node
        # fix(v[b], 1; force=true)        
    end
end



@NLexpression(m, p_ft[l in lineset],
    +v[lines[l].fbus]^2 * yff_r[l]
    +
    v[lines[l].fbus] * v[lines[l].tbus] * (
        +yft_r[l] * cos(θ[lines[l].fbus] - θ[lines[l].tbus])
        +
        yft_i[l] * sin(θ[lines[l].fbus] - θ[lines[l].tbus])
    )
)

@NLexpression(m, q_ft[l in lineset],
    -v[lines[l].fbus]^2 * yff_i[l]
    +
    v[lines[l].fbus] * v[lines[l].tbus] * (
        -yft_i[l] * cos(θ[lines[l].fbus] - θ[lines[l].tbus])
        +
        yft_r[l] * sin(θ[lines[l].fbus] - θ[lines[l].tbus])
    )
)

@NLexpression(m, p_tf[l in lineset],
    +v[lines[l].tbus]^2 * ytt_r[l]
    +
    v[lines[l].tbus] * v[lines[l].fbus] * (
        +ytf_r[l] * cos(θ[lines[l].tbus] - θ[lines[l].fbus])
        +
        ytf_i[l] * sin(θ[lines[l].tbus] - θ[lines[l].fbus])
    )
)

@NLexpression(m, q_tf[l in lineset],
    -v[lines[l].tbus]^2 * ytt_i[l]
    +
    v[lines[l].tbus] * v[lines[l].fbus] * (
        -ytf_i[l] * cos(θ[lines[l].tbus] - θ[lines[l].fbus])
        +
        ytf_r[l] * sin(θ[lines[l].tbus] - θ[lines[l].fbus])
    )
)

@NLconstraint(m, PBalance[b in busset],
    sum(p_ft[l] for l in buses[b].outline)
    + sum(p_tf[l] for l in buses[b].inline)
    +
 - sum(pg[g] for g in B_gn[b])
    + buses[b].Pd
    + buses[b].Gs * v[b]^2
    ==
    0)

@NLconstraint(m, QBalance[b in busset],
    sum(q_ft[l] for l in buses[b].outline)
    + sum(q_tf[l] for l in buses[b].inline)
    +
 - sum(qg[g] for g in B_gn[b])
    +
    buses[b].Qd
    -
    buses[b].Bs * v[b]^2
    ==
    0)

@NLconstraint(m, LineCapFW[l in lineset; lines[l].u != 0], p_ft[l]^2 + q_ft[l]^2 <= (lines[l].u)^2)
@NLconstraint(m, LineCapBW[l in lineset; lines[l].u != 0], p_tf[l]^2 + q_tf[l]^2 <= (lines[l].u)^2)


@constraint(m, AngDiffmax[l in lineset], θ[lines[l].fbus] - θ[lines[l].tbus] <= min(lines[l].angmax, deg2rad(60)))
@constraint(m, AngDiffmin[l in lineset], θ[lines[l].fbus] - θ[lines[l].tbus] >= max(lines[l].angmin, -deg2rad(60)))

@expression(m, OCgen[g in genset],
    generators[g].cost[1] * (pg[g] * datamat["baseMVA"])^2
    + generators[g].cost[2] * (pg[g] * datamat["baseMVA"])
    + generators[g].cost[3])
@expression(m, OC, sum(OCgen[g] for g in genset))

@objective(m, Min, OC)


optimize!(m)
@show termination_status(m)
@show objective_value(m)

