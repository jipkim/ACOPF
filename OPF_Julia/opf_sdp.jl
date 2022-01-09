# SDP formulation is based on `Lavaei, Javad, and Steven H. Low. "Zero duality gap in optimal power flow problem." IEEE Transactions on Power Systems 27, no. 1 (2011): 92-107.`

using JuMP
using Mosek, MosekTools
# using SCS
using DataFrames # for comparison 

include(string(pwd(), "/src/struct_network.jl"))
include(string(pwd(), "/src/readnetwork.jl"))
include(string(pwd(), "/src/ybus.jl"))
include(string(pwd(), "/src/sdp_ymat.jl"))
#----------------------------------------------------------------------------------#
solver = "Mosek"
# solver = "SCS"

testsystem = "case_ieee30"


buses, lines, generators, datamat = readnetwork(testsystem)


nline = length(lines)
nbus = length(buses)
ngen = length(generators)

lineset = 1:nline
busset = 1:nbus
genset = 1:ngen

B_g = []
for g in genset
   push!(B_g, generators[g].location)
end

B_gn = Array{Array{Int64,1}}(undef, nbus)
for ii in 1:nbus
   B_gn[ii] = Int64[]
end
for g in genset
   push!(B_gn[generators[g].location], g)
end

Ybus, yff, yft, ytf, ytt = ybus( buses, lines )
Yk, Yk_, Mk, Ylineft, Ylinetf, Y_lineft, Y_linetf, YL, YL_ = sdp_ymat(lines, Ybus);

if size(datamat["gencost"],2) == 7
    c_k2 =  datamat["gencost"][:,5] * datamat["baseMVA"]^2; 
    c_k1 =  datamat["gencost"][:,6] * datamat["baseMVA"]; 
    c_k0 =  datamat["gencost"][:,7]
elseif size(datamat["gencost"],2) == 6
    c_k2 =  spzeros( size( datamat["gencost"], 1 ), 1 ); 
    c_k1 =  datamat["gencost"][:,5] * datamat["baseMVA"]; 
    c_k0 =  datamat["gencost"][:,6]
end

########################################################################
#---------------------------------------------------------------------------------#
if solver == "SCS"
   m = Model(SCS.Optimizer)   
elseif solver == "Mosek"
   m = Model(Mosek.Optimizer)
end


# p5 right column: a small resistance (1e-5) is added to each transformer ...
# ∴ connected graph
for l in lineset
    if lines[l].r == 0
        # lines[l].r == 1e-5;
    end
end



@variable(m, W[1:2*nbus, 1:2*nbus], PSD) 
@variable(m, α_k[g in genset])

@objective(m, Min, sum(α_k[g] for g in genset))

for b in busset 
    if buses[b].btype == 3        
        @constraint(m, [n in busset], W[n, nbus+b] == 0)
        @constraint(m, [n in busset], W[nbus+b, n] == 0)
    end
end

# Equations (4a)
@constraint(m, [b in busset], sum(generators[g].Pmin for g in B_gn[b]) - buses[b].Pd <= tr( Yk(b) * W ))
@constraint(m, [b in busset], sum(generators[g].Pmax for g in B_gn[b]) - buses[b].Pd >= tr( Yk(b) * W ))

# Equations (4b)
@constraint(m, [b in busset], sum(generators[g].Qmin for g in B_gn[b]) - buses[b].Qd <= tr( Yk_(b) * W ))
@constraint(m, [b in busset], sum(generators[g].Qmax for g in B_gn[b]) - buses[b].Qd >= tr( Yk_(b) * W ))

# Equations (4c)
@constraint(m, V_lb[b in busset], buses[b].Vmin^2 <= tr( Mk(b) * W ))
@constraint(m, V_ub[b in busset], buses[b].Vmax^2 >= tr( Mk(b) * W ))

# Equations (4e)
@constraint(m, LineCap_P_lb[l in lineset; lines[l].u != 0], - lines[l].u <= tr( Ylineft(l) * W ))
@constraint(m, LineCap_P_ub[l in lineset; lines[l].u != 0], lines[l].u >= tr( Ylineft(l) * W ))

# Equations (5): 3x3 matrix
@SDconstraint(m, LineCap_S[l in lineset; lines[l].u != 0], 
    [
        - lines[l].u^2 (tr( Ylineft(l) * W )) (tr( Y_lineft(l) * W ));
        (tr( Ylineft(l) * W )) -1 0;
        (tr( Y_lineft(l) * W )) 0 -1
    ] <= 0
        )

# Equations (6): 2x2 matrix
@SDconstraint(m, [g in genset], 
    [
        c_k1[g] * tr( Yk(B_g[g]) * W ) - α_k[g] + c_k0[g] + c_k1[g] * buses[B_g[g]].Pd     sqrt( c_k2[g] ) * tr( Yk(B_g[g]) * W ) + sqrt( c_k2[g] ) * ( buses[B_g[g]].Pd );
        sqrt( c_k2[g] ) * tr( Yk(B_g[g]) * W ) + sqrt( c_k2[g] ) * ( buses[B_g[g]].Pd )    -1
    ] <= 0
        )



optimize!(m)
@show termination_status(m)
@show objective_value(m)


