clear all; clc;
testsystem = "case_ieee30";
load("data/network/mat/" + testsystem + ".mat"); % load function from MATPOWER

%%
[buses, lines, generators] = readnetwork(mpc);


lineset = 1:size( lines, 1 );
busset = 1:size( buses, 1 );
genset = 1:size( generators, 1);

B_g = [];
for g = genset
    B_g = [ B_g; generators(g).location ];
end

%%
B_gn_cell = cell( size( buses, 1 ), 1 );
for g = genset
    B_gn_cell{generators(g).location} = [B_gn_cell{generators(g).location}; g];
end
B_gn = containers.Map( busset, B_gn_cell);

[ Ybus, yff, yft, ytf, ytt ] = ybus( buses, lines );

yff_r = real(yff);
yff_i = imag(yff);
ytt_r = real(ytt);
ytt_i = imag(ytt);
yft_r = real(yft);
yft_i = imag(yft);
ytf_r = real(ytf);
ytf_i = imag(ytf);

%%
Constraints = [];

theta = sdpvar( length(busset), 1 );
pg = sdpvar( length(genset), 1 );
p_ft = sdpvar( length(lineset), 1 );
p_tf = sdpvar( length(lineset), 1 );
% OC = sdpvar;


for b = busset
    Constraints = [Constraints, - deg2rad(180) <= theta(b) <= deg2rad(180)];
    %     Constraints = [Constraints, - pi <= theta(b) <= pi];
end
for g = genset
    Constraints = [Constraints, generators(g).Pmin <= pg(g) <= generators(g).Pmax];
end

for b = busset
    if buses(b).btype == 3
        Constraints = [Constraints, theta(b) == 0];
    end
end

for l = lineset
    Constraints = [Constraints, p_ft(l) == yft_i(l) * ( theta( lines(l).fbus ) - theta( lines(l).tbus ) )];
    Constraints = [Constraints, p_tf(l) == ytf_i(l) * ( theta( lines(l).tbus ) - theta( lines(l).fbus ) )];
    
    if lines(l).u ~= 0
        Constraints = [Constraints, - lines(l).u <= p_ft(l) <= lines(l).u ];
        Constraints = [Constraints, - lines(l).u <= p_tf(l) <= lines(l).u];
    end
    
    Constraints = [Constraints, lines(l).angmin <= theta(lines(l).fbus) - theta(lines(l).tbus) <= lines(l).angmax];
    % Constraints = [Constraints, max(lines(l).angmin, -deg2rad(60)) <= theta(lines(l).fbus) - theta(lines(l).tbus) <= min(lines(l).angmax, deg2rad(60))];
end


for b = busset
    Constraints = [Constraints, ...
        sum( p_ft( buses(b).outline ) ) ...
        + sum( p_tf( buses(b).inline ) ) ...
        - sum( pg( B_gn(b) )) ...
        + buses(b).Pd ...
        + buses(b).Gs * 1.0^2 ...
        == 0 ...
        ];
end

OC = 0;
for g = genset
    OC = OC + generators(g).cost(1) * (pg(g) * mpc.baseMVA)^2 ...
        + generators(g).cost(2) * (pg(g) * mpc.baseMVA) ...
        + generators(g).cost(3);
end

Objective = OC;





%% Run the optimization
options = sdpsettings('verbose',1,'debug',1,'solver','Gurobi');
% options = sdpsettings('verbose',1,'debug',1,'solver','Mosek');
sol = optimize(Constraints, Objective, options);
% sol = optimize(Constraints, Objective);
yalmip_objvalue = value(OC);
