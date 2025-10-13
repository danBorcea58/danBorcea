import numpy as np

def qck(angle):
    twopi = 2 * np.pi
    diff = twopi * (np.fix(angle / twopi) + min(0, np.sign(angle)))
    return angle - diff

def h_E(E, y, m, Nrev):
    tanE2 = np.tan(E / 2)
    h = (Nrev * np.pi + E - np.sin(E)) / tanE2**3 - 4 / m * (y**3 - y**2)
    dh = (1 - np.cos(E)) / tanE2**3 - 1.5 * (Nrev * np.pi + E - np.sin(E)) * (1 / np.cos(E / 2))**2 / tanE2**4
    return h, dh

def lambertMR_one_revolution(RI, RF, TOF, MU, orbitType):
    TOL = 1e-14
    nitermax = 2000
    TWOPI = 2 * np.pi

    RIM2 = np.dot(RI, RI)
    RIM = np.sqrt(RIM2)
    RFM2 = np.dot(RF, RF)
    RFM = np.sqrt(RFM2)
    CTH = np.dot(RI, RF) / (RIM * RFM)
    CR = np.cross(RI, RF)
    STH = np.linalg.norm(CR) / (RIM * RFM)

    if orbitType == 0 and CR[2] < 0:
        STH = -STH
    elif orbitType == 1 and CR[2] > 0:
        STH = -STH

    THETA = qck(np.arctan2(STH, CTH))
    if THETA == 0 or THETA == TWOPI:
        return 0, 0, 0, 2, np.zeros(3), np.zeros(3), 0, 0

    B1 = np.sign(STH) if STH != 0 else 1
    C = np.sqrt(RIM2 + RFM2 - 2 * RIM * RFM * CTH)
    S = (RIM + RFM + C) / 2
    BETA = 2 * np.arcsin(np.sqrt((S - C) / S))
    PMIN = TWOPI * np.sqrt(S**3 / (8 * MU))
    LAMBDA = B1 * np.sqrt((S - C) / S)

    if 4 * TOF * LAMBDA == 0 or abs((S - C) / S) < TOL:
        return 0, 0, 0, -1, np.zeros(3), np.zeros(3), 0, 0

    if THETA * 180 / np.pi <= 5:
        W = np.arctan((RFM / RIM)**0.25) - np.pi / 4
        R1 = (np.sin(THETA / 4))**2
        S1 = (np.tan(2 * W))**2
        L = (R1 + S1) / (R1 + S1 + np.cos(THETA / 2))
    else:
        L = ((1 - LAMBDA) / (1 + LAMBDA))**2

    M = 8 * MU * TOF**2 / (S**3 * (1 + LAMBDA)**6)
    TPAR = (np.sqrt(2 / MU) / 3) * (S**1.5 - B1 * (S - C)**1.5)
    L1 = (1 - L) / 2

    Y = 1
    X0 = 0 if (TOF - TPAR) <= 1e-3 else L
    X = -1e8
    N = 0

    while abs(X0 - X) >= abs(X) * TOL + TOL and N <= nitermax:
        N += 1
        X = X0
        ETA = X / (np.sqrt(1 + X) + 1)**2

        DELTA = 1
        U = 1
        SIGMA = 1
        M1 = 0
        while abs(U) > TOL and M1 <= nitermax:
            M1 += 1
            GAMMA = (M1 + 3)**2 / (4 * (M1 + 3)**2 - 1)
            DELTA = 1 / (1 + GAMMA * ETA * DELTA)
            U *= (DELTA - 1)
            SIGMA += U

        C1 = 8 * (np.sqrt(1 + X) + 1) / (3 + 1 / (5 + ETA + (9 * ETA / 7) * SIGMA))

        if N == 1:
            DENOM = (1 + 2 * X + L) * (3 * C1 + X * C1 + 4 * X)
            H1 = (L + X)**2 * (C1 + 1 + 3 * X) / DENOM
            H2 = M * (C1 + X - L) / DENOM
        else:
            QR = np.sqrt(L1**2 + M / Y**2)
            XPLL = QR - L1
            LP2XP1 = 2 * QR
            DENOM = LP2XP1 * (3 * C1 + X * C1 + 4 * X)
            H1 = (XPLL**2) * (C1 + 1 + 3 * X) / DENOM
            H2 = M * (C1 + X - L) / DENOM

        B = 27 * H2 / (4 * (1 + H1)**3)
        U = -B / (2 * (np.sqrt(B + 1) + 1))

        DELTA = 1
        U0 = 1
        SIGMA = 1
        N1 = 0
        while N1 < nitermax and abs(U0) >= TOL:
            if N1 == 0:
                GAMMA = 4 / 27
                DELTA = 1 / (1 - GAMMA * U * DELTA)
                U0 *= (DELTA - 1)
                SIGMA += U0
            else:
                for i in range(2):
                    if i == 0:
                        GAMMA = 2 * (3 * N1 + 1) * (6 * N1 - 1) / (9 * (4 * N1 - 1) * (4 * N1 + 1))
                    else:
                        GAMMA = 2 * (3 * N1 + 2) * (6 * N1 + 1) / (9 * (4 * N1 + 1) * (4 * N1 + 3))
                    DELTA = 1 / (1 - GAMMA * U * DELTA)
                    U0 *= (DELTA - 1)
                    SIGMA += U0
            N1 += 1

        KU = (SIGMA / 3)**2
        Y = ((1 + H1) / 3) * (2 + np.sqrt(B + 1) / (1 - 2 * U * KU))
        X0 = np.sqrt(((1 - L) / 2)**2 + M / Y**2) - (1 + L) / 2

    if N1 >= nitermax or N >= nitermax:
        return 0, 0, 0, 4, np.zeros(3), np.zeros(3), 0, 0

    CONST = M * S * (1 + LAMBDA)**2
    A = CONST / (8 * X0 * Y**2)
    R11 = (1 + LAMBDA)**2 / (4 * TOF * LAMBDA)
    S11 = Y * (1 + X0)
    T11 = CONST / S11

    VI = -R11 * (S11 * (RI - RF) - T11 * RI / RIM)
    VF = -R11 * (S11 * (RI - RF) + T11 * RF / RFM)

    P = (2 * RIM * RFM * Y**2 * (1 + X0)**2 * np.sin(THETA / 2)**2) / CONST
    E = np.sqrt(1 - P / A)

    return A, P, E, 0, VI, VF, TPAR, THETA
