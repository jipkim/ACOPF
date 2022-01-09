mutable struct Bus
   bindex::Int
   btype::Int
   Pd::Float64
   Qd::Float64
   Gs::Float64
   Bs::Float64
   area::Int
   Vm::Float64
   Va::Float64
   baseKV::Float64
   bzone::Int
   Vmax::Float64
   Vmin::Float64
   children::Vector{Int}
   ancestor::Vector{Int}
   inline::Vector{Int}
   outline::Vector{Int}
   function Bus(bindex, btype, Pd, Qd, Gs, Bs, area, Vm, Va, baseKV, bzone, Vmax, Vmin)
      b = new(bindex, btype)
      b.Pd = Pd
      b.Qd = Qd
      b.Gs = Gs
      b.Bs = Bs
      b.area = area
      b.Vm = Vm
      b.Va = Va
      b.baseKV = baseKV
      b.bzone = bzone
      b.Vmax = Vmax
      b.Vmin = Vmin
      b.children = Int[]
      b.ancestor = Int[]
      b.inline = Int[]
      b.outline = Int[]
      return b
   end
end

mutable struct Line
   lindex::Int
   fbus::Int # the "from" node
   tbus::Int # the "to" node
   r::Float64 # the resistance value
   x::Float64 # the reactance value
   b::Float64 # the susceptance value
   u::Float64 # the capacity of the line
   shft::Float64 # transformer phase shift (degress)
   tap::Float64
   angmin::Float64
   angmax::Float64 
   function Line(lindex, fbus, tbus, r, x, b, u, angmin, angmax, shft, tap)
      line = new(lindex, fbus, tbus)
      line.r = r
      line.x = x
      line.b = b
      line.u = u
      line.shft = shft
      line.tap = tap
      line.angmin = angmin
      line.angmax = angmax
      return line
   end
end


mutable struct Generator
   gindex::Int
   gtype::String
   location::Int
   Pg::Float64
   Qg::Float64
   Qmax::Float64
   Qmin::Float64
   Vg::Float64
   mBase::Float64
   status::Int
   Pmax::Float64
   Pmin::Float64
   cost::Vector{Float64}
   SUcost::Float64
   SDcost::Float64
   RU::Float64
   RD::Float64
   UPtime::Float64
   DNtime::Float64
   function Generator(gindex, gtype, location, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, cost, SUcost, SDcost, RU, RD, UPtime, DNtime)
      g = new(gindex, gtype, location)
      g.Pg = Pg
      g.Qg = Qg
      g.Qmax = Qmax
      g.Qmin = Qmin
      g.Vg = Vg
      g.mBase = mBase
      g.status = status
      g.Pmax = Pmax
      g.Pmin = Pmin
      g.cost = cost
      g.SUcost = SUcost
      g.SDcost = SDcost
      g.RU = RU
      g.RD = RD
      g.UPtime = UPtime
      g.DNtime = DNtime
      return g
   end
end
