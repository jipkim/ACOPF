function [Ybus, yff, yft, ytf, ytt] = ybus( buses, lines )

nbus = length(buses);
nline = length(lines);

Ybus = zeros(nbus, nbus);
yff = zeros(nline,1);
yft = zeros(nline,1);
ytf = zeros(nline,1);
ytt = zeros(nline,1);
tau = zeros(nline,1);
for l = 1:nline
    if lines(l).tap == 0
        tau(l) = exp( + 1j * lines(l).shft );
    else
        tau(l) = lines(l).tap * exp( + 1j * lines(l).shft );
    end
    
    ytt(l) = inv( lines(l).r + 1j * lines(l).x ) + 1j * lines(l).b / 2 ;
    yff(l) = (inv( lines(l).r + 1j * lines(l).x ) + 1j * lines(l).b / 2 ) / (tau(l) .* conj(tau(l)));
    yft(l) = - ( inv( lines(l).r + 1j * lines(l).x ) ) / conj(tau(l));
    ytf(l) = - ( inv( lines(l).r + 1j * lines(l).x ) ) / tau(l);
    
    Ybus(lines(l).fbus, lines(l).tbus) = yft(l);
    Ybus(lines(l).tbus, lines(l).fbus) = ytf(l);
    Ybus(lines(l).fbus, lines(l).fbus) = Ybus(lines(l).fbus, lines(l).fbus) + yff(l);
    Ybus(lines(l).tbus, lines(l).tbus) = Ybus(lines(l).tbus, lines(l).tbus) + ytt(l);
end
for b = 1:nbus
    Ybus(b,b) = Ybus(b,b) + buses(b).Gs + 1j * buses(b).Bs;
end
Ybus = sparse( Ybus );
end