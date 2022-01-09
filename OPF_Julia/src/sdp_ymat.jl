using SparseArrays # blockdiag() function
using LinearAlgebra # Identity Matrix
function sdp_ymat(lines, Ybus)

nbus = max(maximum([lines[l].fbus for l in lineset]), maximum([lines[l].tbus for l in lineset]))
e = (k) -> sparse(Matrix(I, nbus, nbus)[:,k])
Yk_small = (k) -> e(k)*e(k)'*Ybus

Yk = (k) -> (1/2)*
    [real(Yk_small(k) + conj.(Yk_small(k)')) imag(conj.(Yk_small(k)') - Yk_small(k));
    imag(Yk_small(k) - conj.(Yk_small(k)')) real(Yk_small(k) + conj.(Yk_small(k)'))]    

Yk_ = (k) -> -(1/2)*
    [imag(Yk_small(k) + conj.(Yk_small(k)')) real(Yk_small(k) - conj.(Yk_small(k)'));
    real(conj.(Yk_small(k)') - Yk_small(k)) imag(Yk_small(k) + conj.(Yk_small(k)'))]

Mk = (k) -> blockdiag( e(k)*e(k)', e(k)*e(k)' )

gl = (l) -> real( 1 / (lines[l].r+im*lines[l].x)); # Real part of line admittance
bl = (l) -> imag( 1 / (lines[l].r+im*lines[l].x)); # Imaginary part of line admittance

tau = (l) -> ifelse(lines[l].tap == 0, 1, lines[l].tap); 
theta = (l) -> lines[l].shft
gbcosft = (l) -> gl(l)*cos(theta(l)) + bl(l)*cos(theta(l)+pi/2)
gbsinft = (l) -> gl(l)*sin(theta(l)) + bl(l)*sin(theta(l)+pi/2)
gbcostf = (l) -> gl(l)*cos(-theta(l)) + bl(l)*cos(-theta(l)+pi/2)
gbsintf = (l) -> gl(l)*sin(-theta(l)) + bl(l)*sin(-theta(l)+pi/2)

Ylineft = (l) -> 0.5*(
    sparse(    
        [lines[l].fbus, lines[l].fbus, lines[l].fbus, lines[l].fbus+nbus, lines[l].fbus+nbus, lines[l].fbus+nbus], 
        [lines[l].fbus, lines[l].tbus, lines[l].tbus+nbus, lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus+nbus], 
        [gl(l)/(tau(l)^2), -gbcosft(l)/tau(l), gbsinft(l)/tau(l), gl(l)/(tau(l)^2), -gbsinft(l)/tau(l), -gbcosft(l)/tau(l)],
                2*nbus, 2*nbus) 
    + 
    sparse(    
        [lines[l].fbus, lines[l].fbus, lines[l].fbus, lines[l].fbus+nbus, lines[l].fbus+nbus, lines[l].fbus+nbus],
        [lines[l].fbus, lines[l].tbus, lines[l].tbus+nbus, lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus+nbus],
        [gl(l)/(tau(l)^2), -gbcosft(l)/tau(l), gbsinft(l)/tau(l), gl(l)/(tau(l)^2), -gbsinft(l)/tau(l), -gbcosft(l)/tau(l)],
            2*nbus,2*nbus)')
                                                           
Y_lineft = (l) -> 0.5*(
    sparse(     
        [lines[l].fbus, lines[l].fbus,  lines[l].fbus,  lines[l].fbus+nbus, lines[l].fbus+nbus, lines[l].fbus+nbus ],
        [lines[l].fbus, lines[l].tbus,  lines[l].tbus+nbus, lines[l].fbus+nbus, lines[l].tbus,  lines[l].tbus+nbus ], 
        [-(bl(l)+lines[l].b/2)/(tau(l)^2), gbsinft(l)/tau(l), gbcosft(l)/tau(l), -(bl(l)+lines[l].b/2)/(tau(l)^2), -gbcosft(l)/tau(l), gbsinft(l)/tau(l)],
            2*nbus,2*nbus) 
    + 
    sparse(     
        [lines[l].fbus, lines[l].fbus, lines[l].fbus,  lines[l].fbus+nbus, lines[l].fbus+nbus, lines[l].fbus+nbus ], 
        [lines[l].fbus, lines[l].tbus, lines[l].tbus+nbus, lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus+nbus], 
        [-(bl(l)+lines[l].b/2)/(tau(l)^2), gbsinft(l)/tau(l), gbcosft(l)/tau(l), -(bl(l)+lines[l].b/2)/(tau(l)^2), -gbcosft(l)/tau(l), gbsinft(l)/tau(l)],
            2*nbus,2*nbus)')

Ylinetf = (l) -> 0.5*(
    sparse(    
        [lines[l].fbus, lines[l].fbus,  lines[l].fbus+nbus, lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus+nbus ],
        [lines[l].tbus, lines[l].tbus+nbus, lines[l].tbus,  lines[l].tbus+nbus, lines[l].tbus, lines[l].tbus+nbus ],
        [-gbcostf(l)/tau(l),   -gbsintf(l)/tau(l),    gbsintf(l)/tau(l),  -gbcostf(l)/tau(l), gl(l), gl(l)],
            2*nbus,2*nbus)
    +
    sparse(
        [lines[l].fbus, lines[l].fbus, lines[l].fbus+nbus, lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus+nbus],
        [lines[l].tbus, lines[l].tbus+nbus, lines[l].tbus, lines[l].tbus+nbus, lines[l].tbus, lines[l].tbus+nbus],
        [-gbcostf(l)/tau(l), -gbsintf(l)/tau(l), gbsintf(l)/tau(l), -gbcostf(l)/tau(l), gl(l), gl(l)],
            2*nbus,2*nbus)')
                              
Y_linetf = (l) -> 0.5*(
    sparse(
        [lines[l].fbus, lines[l].fbus, lines[l].fbus+nbus, lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus+nbus],
        [lines[l].tbus, lines[l].tbus+nbus, lines[l].tbus, lines[l].tbus+nbus, lines[l].tbus, lines[l].tbus+nbus],
        [gbsintf(l)/tau(l),   -gbcostf(l)/tau(l),    gbcostf(l)/tau(l),     gbsintf(l)/tau(l),     -(bl(l)+lines[l].b/2), -(bl(l)+lines[l].b/2)],
            2*nbus,2*nbus) 
    +
    sparse(
        [lines[l].fbus, lines[l].fbus, lines[l].fbus+nbus, lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus+nbus],
        [lines[l].tbus, lines[l].tbus+nbus, lines[l].tbus, lines[l].tbus+nbus, lines[l].tbus, lines[l].tbus+nbus],
        [gbsintf(l)/tau(l),   -gbcostf(l)/tau(l),    gbcostf(l)/tau(l),     gbsintf(l)/tau(l),     -(bl(l)+lines[l].b/2),  -(bl(l)+lines[l].b/2)],
                2*nbus,2*nbus)')
                              

YL = (l) -> sparse(
    [lines[l].fbus, lines[l].fbus, lines[l].fbus+nbus, lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus, lines[l].tbus+nbus, lines[l].tbus+nbus],
    [lines[l].fbus, lines[l].tbus, lines[l].fbus+nbus, lines[l].tbus+nbus, lines[l].fbus, lines[l].tbus, lines[l].fbus+nbus, lines[l].tbus+nbus],
    [1, -1, 1, -1, -1, 1, -1, 1],
        2*nbus,2*nbus) * lines[l].r * ( gl(l)^2 + bl(l)^2 )
                              
YL_ = (l) -> (sparse(
    [lines[l].fbus, lines[l].fbus, lines[l].fbus+nbus, lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus, lines[l].tbus+nbus, lines[l].tbus+nbus],
    [lines[l].fbus, lines[l].tbus, lines[l].fbus+nbus, lines[l].tbus+nbus, lines[l].fbus, lines[l].tbus, lines[l].fbus+nbus, lines[l].tbus+nbus],
    [1, -1, 1, -1, -1, 1, -1, 1],
        2*nbus,2*nbus) * lines[l].x * ( gl(l)^2 + bl(l)^2 )
    -
    sparse(
    [lines[l].fbus, lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus+nbus],
    [lines[l].fbus, lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus+nbus],
    [1, 1, 1, 1],
        2*nbus, 2*nbus) * lines[l].b / 2)

    return Yk, Yk_, Mk, Ylineft, Ylinetf, Y_lineft, Y_linetf, YL, YL_
end