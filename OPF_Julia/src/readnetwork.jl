using MAT, DelimitedFiles
function readnetwork(testsystem)
    filename_mat = string(pwd(),"/data/network/mat/$testsystem.mat")
    if isfile(filename_mat)
        mpc = matread(string(pwd(),"/data/network/mat/","$testsystem.mat"))["mpc"]
    end
    ########################################
    datamat = mpc    
    
    if mpc["version"] != "2"
        error("Testsystem data is not compatible. Check matpower data version")
    end
    #########################
    bidxmap = Dict(mpc["bus"][n,1] => n for n = 1:length(mpc["bus"][:,1]))
    buses = Bus[]
    for i in 1:size(mpc["bus"],1)
        bindex = mpc["bus"][i,1]
        btype = mpc["bus"][i,2]
        Pd = mpc["bus"][i,3]/mpc["baseMVA"]
        Qd = mpc["bus"][i,4]/mpc["baseMVA"]
        Gs = mpc["bus"][i,5]/mpc["baseMVA"] # shunt conductance
        Bs = mpc["bus"][i,6]/mpc["baseMVA"] # shunt susceptance
        area = mpc["bus"][i,7]
        Vm = mpc["bus"][i,8]
        Va = deg2rad(mpc["bus"][i,9])
        baseKV = mpc["bus"][i,10]
        bzone = mpc["bus"][i,11]
        Vmax = mpc["bus"][i,12]
        Vmin = mpc["bus"][i,13]
        b = Bus(bindex, btype, Pd, Qd, Gs, Bs, area, Vm, Va, baseKV, bzone, Vmax, Vmin)
        push!(buses, b)
    end

    generators = Generator[]
    for i in 1:size(mpc["gen"],1)
        gindex = i
        gtype = "NotDefined"
        location = bidxmap[mpc["gen"][i,1]]        
        Pg = mpc["gen"][i,2]/mpc["baseMVA"]
        Qg = mpc["gen"][i,3]/mpc["baseMVA"]
        Qmax = mpc["gen"][i,4]/mpc["baseMVA"]
        Qmin = mpc["gen"][i,5]/mpc["baseMVA"]
        Vg = mpc["gen"][i,6]
        mBase = mpc["gen"][i,7]
        status = mpc["gen"][i,8]
        Pmax = mpc["gen"][i,9]/mpc["baseMVA"]
        Pmin = mpc["gen"][i,10]/mpc["baseMVA"]        
        if length(mpc["gencost"][i,5:end]) == 3
            cost = [mpc["gencost"][i,5], mpc["gencost"][i,6], mpc["gencost"][i,7]]
        elseif length(mpc["gencost"][i,5:end]) == 2
            cost = [0, mpc["gencost"][i,5], mpc["gencost"][i,6]]
        else
            error("generator cost format is incompatible")
        end
        SUcost = mpc["gencost"][i,2]
        SDcost = mpc["gencost"][i,3]
        RU = 0
        RD = 0
        UPtime = 0
        DNtime = 0
        g = Generator(gindex, gtype, location, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, cost, SUcost, SDcost, RU, RD, UPtime, DNtime)
        push!(generators, g)
    end

    lines = Line[]
    for i in 1:size(mpc["branch"],1)
        lindex = i
        fbus = bidxmap[Int(mpc["branch"][i,1])]
        tbus = bidxmap[Int(mpc["branch"][i,2])]
        r = mpc["branch"][i,3] # resistance
        x = mpc["branch"][i,4] # reactance
        b = mpc["branch"][i,5] # total line charging susceptance
        u = mpc["branch"][i,6]/mpc["baseMVA"]
        tap = mpc["branch"][i,9]
        shft = deg2rad(mpc["branch"][i,10])
        angmin = deg2rad(mpc["branch"][i,12]) # minimum angle difference
        angmax = deg2rad(mpc["branch"][i,13]) # maximum angle difference
        push!(buses[tbus].inline, lindex)
        push!(buses[fbus].outline, lindex)
        l = Line(lindex, fbus, tbus, r, x, b, u, angmin, angmax, shft, tap)
        push!(lines,l)
    end

    return buses, lines, generators, datamat
end
