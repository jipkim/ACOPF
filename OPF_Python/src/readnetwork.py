import os
from scipy.io import loadmat
import numpy as np
from math import radians

from src.struct_network import Bus, Line, Generator

def readnetwork(testsystem):
    filename_mat = str(os.getcwd() + "/data/network/mat/" + testsystem + ".mat")
    mpc = loadmat(filename_mat)["mpc"]

    linedata = mpc["branch"][0, 0]
    busdata = mpc["bus"][0, 0]
    gendata = mpc["gen"][0, 0]
    gcostdata = mpc["gencost"][0, 0]
    baseMVA = mpc["baseMVA"][0, 0][0, 0]

    nbus = busdata.shape[0]
    nline = linedata.shape[0]
    ngen = gendata.shape[0]

    datamat = mpc


    if mpc["version"] != "2":
        print("Testsystem data is not compatible. Check matpower data version")


    #%%
    buses = []
    for i in np.arange(0, nbus):
        bindex = busdata[i, 1-1]
        btype = busdata[i, 2-1]
        Pd = busdata[i, 3-1]/baseMVA
        Qd = busdata[i, 4-1]/baseMVA
        Gs = busdata[i, 5-1]/baseMVA  # shunt conductance
        Bs = busdata[i, 6-1]/baseMVA  # shunt susceptance
        area = busdata[i, 7-1]
        Vm = busdata[i, 8-1]
        Va = radians(busdata[i, 9-1])
        baseKV = busdata[i, 10-1]
        bzone = busdata[i, 11-1]
        Vmax = busdata[i, 12-1]
        Vmin = busdata[i, 13-1]
        b = Bus(bindex, btype, Pd, Qd, Gs, Bs, area,
                Vm, Va, baseKV, bzone, Vmax, Vmin)
        buses.append(b)

    bidxmap = {buses[i].bindex: i for i in np.arange(0, nbus)}
    #%%
    lines = []
    for i in np.arange(0, nline):
        lindex = i+1
        fbus = bidxmap[ int(linedata[i, 1-1]) ]
        tbus = bidxmap[ int(linedata[i, 2-1]) ]
        r = linedata[i, 3-1]  # resistance
        x = linedata[i, 4-1]  # reactance
        b = linedata[i, 5-1]  # total line charging susceptance
        u = linedata[i, 6-1]/baseMVA
        tap = linedata[i, 9-1]
        shft = radians(linedata[i, 10-1])
        angmin = radians(linedata[i, 12-1])  # minimum angle difference
        angmax = radians(linedata[i, 13-1])  # maximum angle difference
        buses[tbus].inline.append(i)
        buses[fbus].outline.append(i)
        l = Line(lindex, fbus, tbus, r, x, b, u, shft, tap, angmin, angmax)
        lines.append(l)
    #%%
    generators = []
    for i in np.arange(0, ngen):
        gindex = i+1
        gtype = "NotDefined"
        location = bidxmap[gendata[i, 1-1]]
        Pg = gendata[i, 2-1]/baseMVA
        Qg = gendata[i, 3-1]/baseMVA
        Qmax = gendata[i, 4-1]/baseMVA
        Qmin = gendata[i, 5-1]/baseMVA
        Vg = gendata[i, 6-1]
        mBase = gendata[i, 7-1]
        status = gendata[i, 8-1]
        Pmax = gendata[i, 9-1]/baseMVA
        Pmin = gendata[i, 10-1]/baseMVA
        if len(gcostdata[i, 5-1:]) == 3:
            cost = [gcostdata[i, 5-1], gcostdata[i, 6-1], gcostdata[i, 7-1]]
        elif len(gcostdata[i, 5-1:]) == 2:
            cost = [0, gcostdata[i, 5-1], gcostdata[i, 6-1]]
        else:
            print("generator cost format is incompatible")

        SUcost = gcostdata[i, 2-1]
        SDcost = gcostdata[i, 3-1]
        RU = 0
        RD = 0
        UPtime = 0
        DNtime = 0
        g = Generator(gindex, gtype, location, Pg, Qg, Qmax, Qmin, Vg, mBase,
                    status, Pmax, Pmin, cost, SUcost, SDcost, RU, RD, UPtime, DNtime)
        generators.append(g)

    return buses, lines, generators, datamat
