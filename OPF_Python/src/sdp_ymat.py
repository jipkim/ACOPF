#%%
import numpy as np
from scipy import sparse
from scipy.linalg import block_diag

#%%
def sdp_ymat( lines, Ybus ):
    nbus = Ybus.shape[0]    
    nline = len(lines)
    # busset = np.arange(0, nbus)
    # lineset = np.arange(0, nline)
    #%%

    def e(k): return np.eye(nbus)[:, k][np.newaxis]  # size of e(k):  (1, nbus)
    def Yk_small(k): return (e(k).T @ e(k)) @ Ybus

    def Yk(k): return (1/2) * \
        np.block([
            [np.real(Yk_small(k) + Yk_small(k).T), np.imag(Yk_small(k).T - Yk_small(k))],
            [np.imag(Yk_small(k) - Yk_small(k).T), np.real(Yk_small(k) + Yk_small(k).T)]
        ])

    def Yk_(k): return -(1/2) * \
        np.block([
            [np.imag(Yk_small(k) + Yk_small(k).T), np.real(Yk_small(k) - Yk_small(k).T)],
            [np.real(Yk_small(k).T - Yk_small(k)), np.imag(Yk_small(k) + Yk_small(k).T)]
        ])

    def Mk(k): return block_diag(e(k).T @ e(k), e(k).T @ e(k))
    # Real part of line admittance
    def gl(l): return np.real(1 / (lines[l].r+1j*lines[l].x))
    # Imaginary part of line admittance
    def bl(l): return np.imag(1 / (lines[l].r+1j*lines[l].x))

    def tau(l): return 1 if lines[l].tap == 0 else lines[l].tap
    def theta(l): return lines[l].shft
    def gbcosft(l): return gl(l)*np.cos(theta(l)) + bl(l)*np.cos(theta(l)+np.pi/2)
    def gbsinft(l): return gl(l)*np.sin(theta(l)) + bl(l)*np.sin(theta(l)+np.pi/2)
    def gbcostf(l): return gl(l)*np.cos(-theta(l)) + bl(l)*np.cos(-theta(l)+np.pi/2)
    def gbsintf(l): return gl(l)*np.sin(-theta(l)) + bl(l)*np.sin(-theta(l)+np.pi/2)


    #%%
    def Ylineft(l): return 0.5*(
        sparse.coo_matrix((
            [gl(l)/(tau(l)**2), -gbcosft(l)/tau(l), gbsinft(l)/tau(l),
                gl(l)/(tau(l)**2), -gbsinft(l)/tau(l), -gbcosft(l)/tau(l)],
            ([lines[l].fbus, lines[l].fbus, lines[l].fbus, lines[l].fbus +
            nbus, lines[l].fbus+nbus, lines[l].fbus+nbus],
            [lines[l].fbus, lines[l].tbus, lines[l].tbus+nbus,
            lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus+nbus])
        ), shape = (2*nbus, 2*nbus))
        +
        sparse.coo_matrix((
            [gl(l)/(tau(l)**2), -gbcosft(l)/tau(l), gbsinft(l)/tau(l),
                gl(l)/(tau(l)**2), -gbsinft(l)/tau(l), -gbcosft(l)/tau(l)],
            ([lines[l].fbus, lines[l].fbus, lines[l].fbus, lines[l].fbus +
            nbus, lines[l].fbus+nbus, lines[l].fbus+nbus],
            [lines[l].fbus, lines[l].tbus, lines[l].tbus+nbus,
            lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus+nbus])
        ), shape=(2*nbus, 2*nbus)).T
        )



    def Y_lineft(l): return 0.5*(
        sparse.coo_matrix((
            [-(bl(l)+lines[l].b/2)/(tau(l)**2), gbsinft(l)/tau(l), gbcosft(l)/tau(l), -
            (bl(l)+lines[l].b/2)/(tau(l)**2), -gbcosft(l)/tau(l), gbsinft(l)/tau(l)],
            ([lines[l].fbus, lines[l].fbus,  lines[l].fbus,  lines[l].fbus +
                nbus, lines[l].fbus+nbus, lines[l].fbus+nbus],
            [lines[l].fbus, lines[l].tbus,  lines[l].tbus+nbus,
                lines[l].fbus+nbus, lines[l].tbus,  lines[l].tbus+nbus])
        ), shape=(2*nbus, 2*nbus))
        +
        sparse.coo_matrix((
            [-(bl(l)+lines[l].b/2)/(tau(l)**2), gbsinft(l)/tau(l), gbcosft(l)/tau(l), -
            (bl(l)+lines[l].b/2)/(tau(l)**2), -gbcosft(l)/tau(l), gbsinft(l)/tau(l)],
            ([lines[l].fbus, lines[l].fbus, lines[l].fbus,  lines[l].fbus +
                nbus, lines[l].fbus+nbus, lines[l].fbus+nbus],
            [lines[l].fbus, lines[l].tbus, lines[l].tbus+nbus,
                lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus+nbus])
        ), shape=(2*nbus, 2*nbus)).T
        )


    def Ylinetf(l): return 0.5*(
        sparse.coo_matrix((
            [-gbcostf(l)/tau(l),   -gbsintf(l)/tau(l),    gbsintf(l) /
            tau(l),  -gbcostf(l)/tau(l), gl(l), gl(l)],
            ([lines[l].fbus, lines[l].fbus,  lines[l].fbus+nbus,
            lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus+nbus],
            [lines[l].tbus, lines[l].tbus+nbus, lines[l].tbus,
                lines[l].tbus+nbus, lines[l].tbus, lines[l].tbus+nbus])
        ), shape = (2*nbus, 2*nbus))
        +
        sparse.coo_matrix((
            [-gbcostf(l)/tau(l), -gbsintf(l)/tau(l), gbsintf(l) /
            tau(l), -gbcostf(l)/tau(l), gl(l), gl(l)],
            ([lines[l].fbus, lines[l].fbus, lines[l].fbus+nbus,
            lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus+nbus],
            [lines[l].tbus, lines[l].tbus+nbus, lines[l].tbus,
                lines[l].tbus+nbus, lines[l].tbus, lines[l].tbus+nbus])
        ), shape = (2*nbus, 2*nbus)).T
        )

    def Y_linetf(l): return 0.5*(
        sparse.coo_matrix((
            [gbsintf(l)/tau(l),   -gbcostf(l)/tau(l),    gbcostf(l)/tau(l),
            gbsintf(l)/tau(l),     -(bl(l)+lines[l].b/2), -(bl(l)+lines[l].b/2)],
            ([lines[l].fbus, lines[l].fbus, lines[l].fbus+nbus,
            lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus+nbus],
            [lines[l].tbus, lines[l].tbus+nbus, lines[l].tbus,
                lines[l].tbus+nbus, lines[l].tbus, lines[l].tbus+nbus])
        ), shape=(2*nbus, 2*nbus))
        +
        sparse.coo_matrix((
            [gbsintf(l)/tau(l),   -gbcostf(l)/tau(l),    gbcostf(l)/tau(l),
            gbsintf(l)/tau(l),     -(bl(l)+lines[l].b/2),  -(bl(l)+lines[l].b/2)],
            ([lines[l].fbus, lines[l].fbus, lines[l].fbus+nbus,
            lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus+nbus],
            [lines[l].tbus, lines[l].tbus+nbus, lines[l].tbus,
                lines[l].tbus+nbus, lines[l].tbus, lines[l].tbus+nbus])
        ), shape=(2*nbus, 2*nbus)).T
        )



    def YL(l): return sparse.coo_matrix((
            [1, -1, 1, -1, -1, 1, -1, 1],
            ([lines[l].fbus, lines[l].fbus, lines[l].fbus+nbus, lines[l].fbus+nbus,
                lines[l].tbus, lines[l].tbus, lines[l].tbus+nbus, lines[l].tbus+nbus],
            [lines[l].fbus, lines[l].tbus, lines[l].fbus+nbus, lines[l].tbus+nbus,
                lines[l].fbus, lines[l].tbus, lines[l].fbus+nbus, lines[l].tbus+nbus])
        ), shape = (2*nbus, 2*nbus)) * lines[l].r * (gl(l)**2 + bl(l)**2)




    def YL_(l): return (sparse.coo_matrix((
            [1, -1, 1, -1, -1, 1, -1, 1],
            ([lines[l].fbus, lines[l].fbus, lines[l].fbus+nbus, lines[l].fbus+nbus,
                lines[l].tbus, lines[l].tbus, lines[l].tbus+nbus, lines[l].tbus+nbus],
            [lines[l].fbus, lines[l].tbus, lines[l].fbus+nbus, lines[l].tbus+nbus,
                lines[l].fbus, lines[l].tbus, lines[l].fbus+nbus, lines[l].tbus+nbus])
        ), shape = (2*nbus, 2*nbus)) * lines[l].x * (gl(l)**2 + bl(l)**2)
        -
        sparse.coo_matrix((
            [1, 1, 1, 1],
            ([lines[l].fbus, lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus+nbus],
            [lines[l].fbus, lines[l].fbus+nbus, lines[l].tbus, lines[l].tbus+nbus])
        ), shape = (2*nbus, 2*nbus)) * lines[l].b / 2)

    return Yk, Yk_, Mk, Ylineft, Ylinetf, Y_lineft, Y_linetf, YL, YL_

