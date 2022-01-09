function [buses, lines, generators] = readnetwork(mpc)
%%
%
if mpc.version ~= "2"
    error("Testsystem data is not compatible. Check matpower data version")
end
bidxmap = containers.Map( mpc.bus(:,1), 1:length( mpc.bus(:,1) ) );

%%
buses = [];
for i = 1:size( mpc.bus, 1 )
    b.bindex = mpc.bus(i,1);
    b.btype = mpc.bus(i,2);
    b.Pd = mpc.bus(i,3)/mpc.baseMVA;
    b.Qd = mpc.bus(i,4)/mpc.baseMVA;
    b.Gs = mpc.bus(i,5)/mpc.baseMVA; % shunt conductance
    b.Bs = mpc.bus(i,6)/mpc.baseMVA; % shunt susceptance
    b.area = mpc.bus(i,7);
    b.Vm = mpc.bus(i,8);
    b.Va = deg2rad(mpc.bus(i,9));
    b.baseKV = mpc.bus(i,10);
    b.bzone = mpc.bus(i,11);
    b.Vmax = mpc.bus(i,12);
    b.Vmin = mpc.bus(i,13);
    b.inline = [];
    b.outline = [];
    buses = [buses; b];
end
%%
generators = [];
for i = 1:size( mpc.gen, 1 )
    g.gindex = i;
    g.gtype = "NotDefined";
    g.location = bidxmap(mpc.gen(i,1));
    g.Pg = mpc.gen(i,2)/mpc.baseMVA;
    g.Qg = mpc.gen(i,3)/mpc.baseMVA;
    g.Qmax = mpc.gen(i,4)/mpc.baseMVA;
    g.Qmin = mpc.gen(i,5)/mpc.baseMVA;
    g.Vg = mpc.gen(i,6);
    g.mBase = mpc.gen(i,7);
    g.status = mpc.gen(i,8);
    g.Pmax = mpc.gen(i,9)/mpc.baseMVA;
    g.Pmin = mpc.gen(i,10)/mpc.baseMVA;
    if length(mpc.gencost(i,5:end)) == 3
        g.cost = [mpc.gencost(i,5), mpc.gencost(i,6), mpc.gencost(i,7)];
    elseif length(mpc.gencost(i,5:end)) == 2
        g.cost = [0, mpc.gencost(i,5), mpc.gencost(i,6)];
    else
        error("generator cost format is incompatible")
    end
    g.SUcost = mpc.gencost(i,2);
    g.SDcost = mpc.gencost(i,3);
    g.RU = 0;
    g.RD = 0;
    g.UPtime = 0;
    g.DNtime = 0;
    generators = [generators; g];
end


%%
lines = [];
for i = 1:size(mpc.branch,1)
    l.lindex = i;
    l.fbus = bidxmap(mpc.branch(i,1));
    l.tbus = bidxmap(mpc.branch(i,2));
    l.r = mpc.branch(i,3); % resistance
    l.x = mpc.branch(i,4); % reactance
    l.b = mpc.branch(i,5); % total line charging susceptance
    l.u = mpc.branch(i,6)/mpc.baseMVA;
    l.tap = mpc.branch(i,9);
    l.shft = deg2rad(mpc.branch(i,10));
    l.angmin = deg2rad(mpc.branch(i,12)); % minimum angle difference
    l.angmax = deg2rad(mpc.branch(i,13)); % maximum angle difference
    buses(l.tbus).inline = [buses(l.tbus).inline, l.lindex];
    buses(l.fbus).outline = [buses(l.fbus).outline, l.lindex];
    lines = [lines; l];
end