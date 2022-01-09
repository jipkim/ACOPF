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
gencost = zeros(length(genset),3);
for g = genset
    gencost(g,:) = generators(g).cost;
end

cvx_solver Gurobi
% cvx_solver Knitro
% cvx_solver Mosek

% Create and solve problem
cvx_begin
variable theta( length(busset), 1 );
variable pg( length(genset), 1 );
variable p_ft( length(lineset), 1 );
variable p_tf( length(lineset), 1 );

OCgen = ( gencost(:,1) .* (pg * mpc.baseMVA).^2 ...
    + gencost(:,2) .* (pg * mpc.baseMVA) ...
    + gencost(:,3) );
OC = sum( OCgen );

minimize( OC )
subject to
for b = busset
    - deg2rad(180) <= theta(b) <= deg2rad(180)
end
for g = genset
    generators(g).Pmin <= pg(g) <= generators(g).Pmax;
end

for b = busset
    if buses(b).btype == 3
        theta(b) == 0;
    end
end
for l = lineset
    p_ft(l) == yft_i(l) * ( theta( lines(l).fbus ) - theta( lines(l).tbus ) );
    p_tf(l) == ytf_i(l) * ( theta( lines(l).tbus ) - theta( lines(l).fbus ) );
    
    if lines(l).u ~= 0
        - lines(l).u <= p_ft(l) <= lines(l).u;
        - lines(l).u <= p_tf(l) <= lines(l).u;
    end
    
    lines(l).angmin <= theta(lines(l).fbus) - theta(lines(l).tbus) <= lines(l).angmax;
    % max(lines(l).angmin, -deg2rad(60)) <= theta(lines(l).fbus) - theta(lines(l).tbus) <= min(lines(l).angmax, deg2rad(60));
end
for b = busset
    sum( p_ft( buses(b).outline ) ) ...
        + sum( p_tf( buses(b).inline ) ) ...
        - sum( pg( B_gn(b) )) ...
        + buses(b).Pd ...
        + buses(b).Gs * 1.0^2 ...
        == 0;
end
cvx_end

%% Run the optimization
cvx_objvalue = value( OC );
%%
