clear all
testsystem = "case_ieee30";
load("data/network/mat/" + testsystem + ".mat"); % load function from MATPOWER




[buses, lines, generators] = readnetwork(mpc);


nline = length(lines);
nbus = length(buses);
ngen = length(generators);

lineset = 1:nline;
busset = 1:nbus;
genset = 1:ngen;


B_g = [];
for g = genset
    B_g = [ B_g; generators(g).location ];
end

%%
B_gn_cell = cell( nbus, 1 );
for g = genset
    B_gn_cell{generators(g).location} = [B_gn_cell{generators(g).location}; g];
end
B_gn = containers.Map( busset, B_gn_cell);


[ Ybus, yff, yft, ytf, ytt ] = ybus( buses, lines );


% Add resistance to transformers/lines -> connected resistive graph
% for l = lineset
%     if lines(l).r == 0
%         lines(l).r = 0.0001;
%     end
% end

[Yk,Yk_,Mk,Ylineft,Ylinetf,Y_lineft,Y_linetf,YL,YL_] = sdp_ymat( lines, Ybus );



cvx_solver Mosek
cvx_begin SDP
variable W( 2*nbus, 2*nbus ) semidefinite
variable alpha_k( ngen, 1 )
minimize sum( alpha_k, 1 )
% slack bus
for b = busset
    if buses(b).btype == 3
        W( :, nbus+b ) == 0;
        W( nbus+b, : ) == 0;
        slack_bus = b;
        display("slack bus: " + slack_bus)
    end
end

for b = busset
    sum([generators(B_gn(b)).Pmin]) - buses(b).Pd <= trace( Yk(b) * W ) <= sum([generators(B_gn(b)).Pmax]) - buses(b).Pd
    sum([generators(B_gn(b)).Qmin]) - buses(b).Qd <= trace( Yk_(b) * W ) <= sum([generators(B_gn(b)).Qmax]) - buses(b).Qd
end

% Equation (4c)
for b = busset
    ( buses(b).Vmin )^2 <= trace( Mk(b) * W ) <= ( buses(b).Vmax )^2
end
% Equation (4e)
for l = lineset
    if lines(l).u ~= 0
        - lines(l).u <= trace( Ylineft(l) * W ) <= lines(l).u
    end
end
% Equations (5)
for l = lineset
    if lines(l).u ~= 0
        [...
            -( lines(l).u )^2, trace( Ylineft(l) * W ), trace( Y_lineft(l) * W );...
            trace( Ylineft(l) * W ), -1, 0;...
            trace( Y_lineft(l) * W ), 0, -1] <= 0
    end
end
% Equation (6)
for g = genset
    [...
        generators(g).cost(2) * ( trace( Yk(B_g(g)) * W ) * mpc.baseMVA ) - alpha_k( g, 1 ) + generators(g).cost(3) + generators(g).cost(2) * ( buses(B_g(g)).Pd * mpc.baseMVA ), ...
        sqrt( generators(g).cost(1) ) * ( trace( Yk(B_g(g)) * W ) * mpc.baseMVA ) + sqrt( generators(g).cost(1) ) * ( buses(B_g(g)).Pd * mpc.baseMVA ); ...
        sqrt( generators(g).cost(1) ) * ( trace( Yk(B_g(g)) * W ) * mpc.baseMVA ) + sqrt( generators(g).cost(1) ) * ( buses(B_g(g)).Pd * mpc.baseMVA ), -1] <= 0
end
cvx_end

cvx_objvalue = sum(alpha_k,1);

%%
