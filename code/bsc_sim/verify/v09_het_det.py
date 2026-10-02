import numpy as np, vlib as V
B, G = V.BETA, V.GAMMA
def twozone(c, D0=1.0, N=100, variant='zoneD'):
    nx, ny = 20, 13
    acc = np.ones((ny, nx), bool)
    isC = np.zeros((ny, nx), bool); isC[:, 7:13] = True
    q = np.where(isC, 1.2 * (1 - c), 1.2 * (1 + c * 0.3 / 0.7))
    D = np.where(isC, 2.0, 0.3) * D0
    if variant == 'uniformD': D = np.full((ny, nx), 0.81 * D0)
    return V.Room(acc, D, q, N)
if __name__ == '__main__':
    ct = 1 - 0.5 / 1.2
    for D0 in (0.01, 0.1, 1, 10, 100, 1000):
        r = twozone(ct, D0)
        print('D0', D0, 'R0/(beta/gamma) = %.4f' % (V.rho(V.ngm(r)) / (B / G)), ' qmean', r.q.mean(), 'qW,qC', r.q.max(), r.q.min())
    for c in (0, 0.2, 0.4, ct, 0.8, 1.0):
        print('c=%.3f R0(D0=1)=%.3f R0(D0=10)=%.3f  uniformD: %.3f' % (c, V.rho(V.ngm(twozone(c, 1))), V.rho(V.ngm(twozone(c, 10))), V.rho(V.ngm(twozone(c, 1, variant='uniformD')))))
    for c in (-0.5, -1, -2):
        print('c=%.1f hot corridor R0(D0=1)=%.3f' % (c, V.rho(V.ngm(twozone(c, 1)))))
