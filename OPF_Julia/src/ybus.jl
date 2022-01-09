using SparseArrays # spzeros
function ybus( buses, lines )

nbus = length(buses)
nline = length(lines)
busset = 1:nbus
lineset = 1:nline

Ybus = spzeros(ComplexF64, nbus, nbus)
yff = spzeros(ComplexF64, nline, 1)
yft = spzeros(ComplexF64, nline, 1)
ytf = spzeros(ComplexF64, nline, 1)
ytt = spzeros(ComplexF64, nline, 1)
tau = spzeros(ComplexF64, nline, 1)
for l in lineset
    tau[l] = ifelse( lines[l].tap == 0, exp( + im * lines[l].shft), lines[l].tap * exp( + im * lines[l].shft) )

    ytt[l] = inv(lines[l].r + im * lines[l].x) + im * lines[l].b / 2 
    yff[l] = (inv(lines[l].r + im * lines[l].x) + im * lines[l].b / 2 ) / (tau[l] .* conj(tau[l]))
    yft[l] = - (inv(lines[l].r + im * lines[l].x)) / conj(tau[l])
    ytf[l] = - (inv(lines[l].r + im * lines[l].x)) / tau[l]

    Ybus[lines[l].fbus, lines[l].tbus] = yft[l]
    Ybus[lines[l].tbus, lines[l].fbus] = ytf[l]
    Ybus[lines[l].fbus, lines[l].fbus] += yff[l]
    Ybus[lines[l].tbus, lines[l].tbus] += ytt[l]
end
for b in busset
    Ybus[b,b] += buses[b].Gs + im * buses[b].Bs    
end

return Ybus, yff, yft, ytf, ytt
end