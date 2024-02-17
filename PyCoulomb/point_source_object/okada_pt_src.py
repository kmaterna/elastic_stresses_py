import numpy as np

# Common blocks
ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD, S2D, C2D = [0.0] * 12
P, Q, S, T, XY, X2, Y2, D2, R, R2, R3, R5, R7, A3, A5, B3, C3, QR, QRX, UY, VY, WY, UZ, VZ, WZ = [0.0] * 25

def initialize_globals():
    global ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD, S2D, C2D
    global P, Q, S, T, XY, X2, Y2, D2, R, R2, R3, R5, R7, A3, A5, B3, C3, QR, QRX, UY, VY, WY, UZ, VZ, WZ

    ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD, S2D, C2D = [0.0] * 12
    P, Q, S, T, XY, X2, Y2, D2, R, R2, R3, R5, R7, A3, A5, B3, C3, QR, QRX, UY, VY, WY, UZ, VZ, WZ = [0.0] * 25
    return


def UA0(X, Y, D, POT1, POT2, POT3, POT4):
    # Function to calculate displacement and strain at depth due to a buried point source (Part-A)

    # Constants
    F0, F1, F3 = 0.0, 1.0, 3.0
    PI2 = 6.283185307179586

    # Initialize Arrays
    U = np.zeros(12)
    DU = np.zeros(12)

    # Strike-slip contribution
    if POT1 != F0:
        DU[0] = ALP1 * Q / R ** 3 + ALP2 * X ** 2 * QR
        DU[1] = ALP1 * X / R ** 3 * SD + ALP2 * XY * QR
        DU[2] = -ALP1 * X / R ** 3 * CD + ALP2 * X * D * QR
        DU[3] = X * QR * (-ALP1 + ALP2 * (F1 + A5))
        DU[4] = ALP1 * A3 / R ** 3 * SD + ALP2 * Y * QR * A5
        DU[5] = -ALP1 * A3 / R ** 3 * CD + ALP2 * D * QR * A5
        DU[6] = ALP1 * (SD / R ** 3 - Y * QR) + ALP2 * F3 * X ** 2 / R ** 5 * UY
        DU[7] = F3 * X / R ** 5 * (-ALP1 * Y * SD + ALP2 * (Y * UY + Q))
        DU[8] = F3 * X / R ** 5 * (ALP1 * Y * CD + ALP2 * D * UY)
        DU[9] = ALP1 * (CD / R ** 3 + D * QR) + ALP2 * F3 * X ** 2 / R ** 5 * UZ
        DU[10] = F3 * X / R ** 5 * (ALP1 * D * SD + ALP2 * Y * UZ)
        DU[11] = F3 * X / R ** 5 * (-ALP1 * D * CD + ALP2 * (D * UZ - Q))
        U += POT1 / PI2 * DU

    # Dip-slip contribution
    if POT2 != F0:
        DU[0] = ALP2 * X * P * QR
        DU[1] = ALP1 * S / R ** 3 + ALP2 * Y * P * QR
        DU[2] = -ALP1 * T / R ** 3 + ALP2 * D * P * QR
        DU[3] = ALP2 * P * QR * A5
        DU[4] = -ALP1 * F3 * X * S / R ** 5 - ALP2 * Y * P * QRX
        DU[5] = ALP1 * F3 * X * T / R ** 5 - ALP2 * D * P * QRX
        DU[6] = ALP2 * F3 * X / R ** 5 * VY
        DU[7] = ALP1 * (S2D / R ** 3 - F3 * Y * S / R ** 5) + ALP2 * (F3 * Y / R ** 5 * VY + P * QR)
        DU[8] = -ALP1 * (C2D / R ** 3 - F3 * Y * T / R ** 5) + ALP2 * F3 * D / R ** 5 * VY
        DU[9] = ALP2 * F3 * X / R ** 5 * VZ
        DU[10] = ALP1 * (C2D / R ** 3 + F3 * D * S / R ** 5) + ALP2 * F3 * Y / R ** 5 * VZ
        DU[11] = ALP1 * (S2D / R ** 3 - F3 * D * T / R ** 5) + ALP2 * (F3 * D / R ** 5 * VZ - P * QR)
        U += POT2 / PI2 * DU

    # Tensile-fault contribution
    if POT3 != F0:
        DU[0] = ALP1 * X / R ** 3 - ALP2 * X * Q * QR
        DU[1] = ALP1 * T / R ** 3 - ALP2 * Y * Q * QR
        DU[2] = ALP1 * S / R ** 3 - ALP2 * D * Q * QR
        DU[3] = ALP1 * A3 / R ** 3 - ALP2 * Q * QR * A5
        DU[4] = -ALP1 * F3 * X * T / R ** 5 + ALP2 * Y * Q * QRX
        DU[5] = -ALP1 * F3 * X * S / R ** 5 + ALP2 * D * Q * QRX
        DU[6] = -ALP1 * F3 * XY / R ** 5 - ALP2 * X * QR * WY
        DU[7] = ALP1 * (C2D / R ** 3 - F3 * Y * T / R ** 5) - ALP2 * (Y * WY + Q) * QR
        DU[8] = ALP1 * (S2D / R ** 3 - F3 * Y * S / R ** 5) - ALP2 * D * QR * WY
        DU[9] = ALP1 * F3 * X * D / R ** 5 - ALP2 * X * QR * WZ
        DU[10] = -ALP1 * (S2D / R ** 3 - F3 * D * T / R ** 5) - ALP2 * Y * QR * WZ
        DU[11] = ALP1 * (C2D / R ** 3 + F3 * D * S / R ** 5) - ALP2 * (D * WZ - Q) * QR
        U += POT3 / PI2 * DU

    # Inflate source contribution
    if POT4 != F0:
        DU[0:3] = -ALP1 * np.array([X, Y, D]) / R ** 3
        DU[3] = -ALP1 * A3 / R ** 3
        DU[4:6] = ALP1 * F3 * np.array([XY, X * D]) / R ** 5
        DU[6:9] = DU[4:7]
        DU[9:12] = ALP1 * np.array([C3, F3 * D * S, -C3]) / R ** 3
        U += POT4 / PI2 * DU

    return U


def UB0(X, Y, D, Z, POT1, POT2, POT3, POT4):
    # Function to calculate displacement and strain at depth due to a buried point source (Part-B)

    # Constants
    F0, F1, F2, F3, F4, F5, F8, F9 = 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 8.0, 9.0
    PI2 = 6.283185307179586

    # Initialize Arrays
    U = np.zeros(12)
    DU = np.zeros(12)

    # Calculate common terms
    C = D + Z
    RD = R + D
    D12 = F1 / (R * RD * RD)
    D32 = D12 * (F2 * R + D) / R2
    D33 = D12 * (F3 * R + D) / (R2 * RD)
    D53 = D12 * (F8 * R2 + F9 * R * D + F3 * D2) / (R2 * R2 * RD)
    D54 = D12 * (F5 * R2 + F4 * R * D + D2) / R3 * D12

    FI1 = Y * (D12 - X2 * D33)
    FI2 = X * (D12 - Y2 * D33)
    FI3 = X / R3 - FI2
    FI4 = -XY * D32
    FI5 = F1 / (R * RD) - X2 * D32
    FJ1 = -F3 * XY * (D33 - X2 * D54)
    FJ2 = F1 / R3 - F3 * D12 + F3 * X2 * Y2 * D54
    FJ3 = A3 / R3 - FJ2
    FJ4 = -F3 * XY / R5 - FJ1
    FK1 = -Y * (D32 - X2 * D53)
    FK2 = -X * (D32 - Y2 * D53)
    FK3 = -F3 * X * D / R5 - FK2

    # Strike-slip contribution
    if POT1 != F0:
        DU[0] = -X2 * QR - ALP3 * FI1 * SD
        DU[1] = -XY * QR - ALP3 * FI2 * SD
        DU[2] = -C * X * QR - ALP3 * FI4 * SD
        DU[3] = -X * QR * (F1 + A5) - ALP3 * FJ1 * SD
        DU[4] = -Y * QR * A5 - ALP3 * FJ2 * SD
        DU[5] = -C * QR * A5 - ALP3 * FK1 * SD
        DU[6] = -F3 * X2 / R5 * UY - ALP3 * FJ2 * SD
        DU[7] = -F3 * XY / R5 * UY - X * QR - ALP3 * FJ4 * SD
        DU[8] = -F3 * C * X / R5 * UY - ALP3 * FK2 * SD
        DU[9] = -F3 * X2 / R5 * UZ + ALP3 * FK1 * SD
        DU[10] = -F3 * XY / R5 * UZ + ALP3 * FK2 * SD
        DU[11] = F3 * X / R5 * (-C * UZ + ALP3 * Y * SD)
        U += POT1 / PI2 * DU

    # Dip-slip contribution
    if POT2 != F0:
        DU[0] = -X * P * QR + ALP3 * FI3 * SDCD
        DU[1] = -Y * P * QR + ALP3 * FI1 * SDCD
        DU[2] = -C * P * QR + ALP3 * FI5 * SDCD
        DU[3] = -P * QR * A5 + ALP3 * FJ3 * SDCD
        DU[4] = Y * P * QRX + ALP3 * FJ1 * SDCD
        DU[5] = C * P * QRX + ALP3 * FK3 * SDCD
        DU[6] = -F3 * X / R5 * VY + ALP3 * FJ1 * SDCD
        DU[7] = -F3 * Y / R5 * VY - P * QR + ALP3 * FJ2 * SDCD
        DU[8] = -F3 * C / R5 * VY + ALP3 * FK1 * SDCD
        DU[9] = -F3 * X / R5 * VZ - ALP3 * FK3 * SDCD
        DU[10] = -F3 * Y / R5 * VZ - ALP3 * FK1 * SDCD
        DU[11] = -F3 * C / R5 * VZ + ALP3 * A3 / R3 * SDCD
        U += POT2 / PI2 * DU

    # Tensile-fault contribution
    if POT3 != F0:
        DU[0] = X * Q * QR - ALP3 * FI3 * SDSD
        DU[1] = Y * Q * QR - ALP3 * FI1 * SDSD
        DU[2] = C * Q * QR - ALP3 * FI5 * SDSD
        DU[3] = Q * QR * A5 - ALP3 * FJ3 * SDSD
        DU[4] = -Y * Q * QRX - ALP3 * FJ1 * SDSD
        DU[5] = -C * Q * QRX - ALP3 * FK3 * SDSD
        DU[6] = X * QR * WY - ALP3 * FJ1 * SDSD
        DU[7] = QR * (Y * WY + Q) - ALP3 * FJ2 * SDSD
        DU[8] = C * QR * WY - ALP3 * FK1 * SDSD
        DU[9] = X * QR * WZ + ALP3 * FK3 * SDSD
        DU[10] = Y * QR * WZ + ALP3 * FK1 * SDSD
        DU[11] = C * QR * WZ - ALP3 * A3 / R3 * SDSD
        U += POT3 / PI2 * DU

    # Inflate source contribution
    if POT4 != F0:
        DU[0:3] = ALP3 * np.array([X, Y, D]) / R ** 3
        DU[3] = ALP3 * A3 / R ** 3
        DU[4:6] = -ALP3 * F3 * np.array([XY, X * D]) / R ** 5
        DU[6:9] = DU[4:7]
        DU[9:12] = -ALP3 * np.array([C3, F3 * D * S, -C3]) / R ** 3
        U += POT4 / PI2 * DU

    return U


def UC0(X, Y, D, Z, POT1, POT2, POT3, POT4):
    # Function to calculate displacement and strain at depth due to a buried point source (Part-B)

    # Constants
    F0, F1, F2, F3, F5, F7, F10, F15 = 0.0, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0
    PI2 = 6.283185307179586

    # Arrays
    U = np.zeros(12)
    DU = np.zeros(12)

    # Initialize variables
    C = D + Z
    Q2 = Q * Q
    R7 = R5 * R2
    A7 = F1 - F7 * X2 / R2
    B5 = F1 - F5 * Y2 / R2
    B7 = F1 - F7 * Y2 / R2
    C5 = F1 - F5 * D2 / R2
    C7 = F1 - F7 * D2 / R2
    D7 = F2 - F7 * Q2 / R2
    QR5 = F5 * Q / R2
    QR7 = F7 * Q / R2
    DR5 = F5 * D / R2

    # Strike-slip contribution
    if POT1 != F0:
        DU[0] = -ALP4 * A3 / R3 * CD + ALP5 * C * QR * A5
        DU[1] = F3 * X / R5 * (ALP4 * Y * CD + ALP5 * C * (SD - Y * QR5))
        DU[2] = F3 * X / R5 * (-ALP4 * Y * SD + ALP5 * C * (CD + D * QR5))
        DU[3] = ALP4 * F3 * X / R5 * (F2 + A5) * CD - ALP5 * C * QRX * (F2 + A7)
        DU[4] = F3 / R5 * (ALP4 * Y * A5 * CD + ALP5 * C * (A5 * SD - Y * QR5 * A7))
        DU[5] = F3 / R5 * (-ALP4 * Y * A5 * SD + ALP5 * C * (A5 * CD + D * QR5 * A7))
        DU[6:8] = DU[4:6]
        DU[8] = F3 * X / R5 * (ALP4 * B5 * CD - ALP5 * F5 * C / R2 * (F2 * Y * SD + Q * B7))
        DU[9] = F3 * X / R5 * (-ALP4 * B5 * SD + ALP5 * F5 * C / R2 * (D * B7 * SD - Y * C7 * CD))
        DU[10] = F3 / R5 * (-ALP4 * D * A5 * CD + ALP5 * C * (A5 * CD + D * QR5 * A7))
        DU[11] = F15 * X / R7 * (ALP4 * Y * D * CD + ALP5 * C * (D * B7 * SD - Y * C7 * CD))
        DU[12] = F15 * X / R7 * (-ALP4 * Y * D * SD + ALP5 * C * (F2 * D * CD - Q * C7))
        U += POT1 / PI2 * DU

    # Dip-slip contribution
    if POT2 != F0:
        DU[0] = ALP4 * F3 * X * T / R5 - ALP5 * C * P * QRX
        DU[1] = -ALP4 / R3 * (C2D - F3 * Y * T / R2) + ALP5 * F3 * C / R5 * (S - Y * P * QR5)
        DU[2] = -ALP4 * A3 / R3 * SDCD + ALP5 * F3 * C / R5 * (T + D * P * QR5)
        DU[3] = ALP4 * F3 * T / R5 * A5 - ALP5 * F5 * C * P * QR / R2 * A7
        DU[4] = F3 * X / R5 * (ALP4 * (C2D - F5 * Y * T / R2) - ALP5 * F5 * C / R2 * (S - Y * P * QR7))
        DU[5] = F3 * X / R5 * (ALP4 * (F2 + A5) * SDCD - ALP5 * F5 * C / R2 * (T + D * P * QR7))
        DU[6:8] = DU[4:6]
        DU[8] = F3 / R5 * (ALP4 * (F2 * Y * C2D + T * B5)
                           + ALP5 * C * (S2D - F10 * Y * S / R2 - P * QR5 * B7))
        DU[9] = F3 / R5 * (ALP4 * Y * A5 * SDCD - ALP5 * C * ((F3 + A5) * C2D + Y * P * DR5 * QR7))
        DU[10] = F3 * X / R5 * (-ALP4 * (S2D - T * DR5) - ALP5 * F5 * C / R2 * (T + D * P * QR7))
        DU[11] = F3 / R5 * (-ALP4 * (D * B5 * C2D + Y * C5 * S2D)
                            - ALP5 * C * ((F3 + A5) * C2D + Y * P * DR5 * QR7))
        DU[12] = F3 / R5 * (-ALP4 * D * A5 * SDCD - ALP5 * C * (S2D - F10 * D * T / R2 + P * QR5 * C7))
        U += POT2 / PI2 * DU

    # Tensile-fault contribution
    if POT3 != F0:
        DU[0] = F3 * X / R5 * (-ALP4 * S + ALP5 * (C * Q * QR5 - Z))
        DU[1] = ALP4 / R3 * (S2D - F3 * Y * S / R2) + ALP5 * F3 / R5 * (
                C * (T - Y + Y * Q * QR5) - Y * Z)
        DU[2] = -ALP4 / R3 * (F1 - A3 * SDSD) - ALP5 * F3 / R5 * (C * (S - D + D * Q * QR5) - D * Z)
        DU[3] = -ALP4 * F3 * S / R5 * A5 + ALP5 * (C * QR * QR5 * A7 - F3 * Z / R5 * A5)
        DU[4] = F3 * X / R5 * (-ALP4 * (S2D - F5 * Y * S / R2)
                               - ALP5 * F5 / R2 * (C * (T - Y + Y * Q * QR7) - Y * Z))
        DU[5] = F3 * X / R5 * (
                ALP4 * (F1 - (F2 + A5) * SDSD) + ALP5 * F5 / R2 * (C * (S - D + D * Q * QR7) - D * Z))
        DU[6:8] = DU[4:6]
        DU[8] = F3 / R5 * (-ALP4 * (F2 * Y * S2D + S * B5)
                           - ALP5 * (C * (F2 * SDSD + F10 * Y * (T - Y) / R2 - Q * QR5 * B7) + Z * B5))
        DU[9] = F3 / R5 * (
                ALP4 * Y * (F1 - A5 * SDSD) + ALP5 * (C * (F3 + A5) * S2D - Y * DR5 * (C * D7 + Z)))
        DU[10] = F3 * X / R5 * (-ALP4 * (C2D + S * DR5)
                                + ALP5 * (F5 * C / R2 * (S - D + D * Q * QR7) - F1 - Z * DR5))
        DU[11] = F3 / R5 * (
                -ALP4 * (D * B5 * S2D - Y * C5 * C2D) - ALP5 * (
                  C * ((F3 + A5) * S2D - Y * DR5 * D7) - Y * (F1 + Z * DR5)))
        DU[12] = F3 / R5 * (
                -ALP4 * D * (F1 - A5 * SDSD) - ALP5 * (
                  C * (C2D + F10 * D * (S - D) / R2 - Q * QR5 * C7) + Z * (F1 + C5)))
        U += POT3 / PI2 * DU

    # Inflate source contribution
    if POT4 != F0:
        DU[0:3] = ALP4 * np.array([F3 * X / R5 * D, F3 * Y / R5 * D, C3 / R3])
        DU[3] = ALP4 * F3 * D / R5 * A5
        DU[4:6] = -ALP4 * F15 * np.array([XY, F3 * X / R5 * C5]) * D / R7
        DU[6:9] = DU[4:7]
        DU[9:12] = ALP4 * F3 * D / R5 * np.array([F2 + C5, -F3 * Y / R5 * C5, F3 * D / R5 * (F2 + C5)])
        U += POT4 / PI2 * DU

    return U


def DCCON0(ALPHA, DIP):
    # Function to calculate medium constants and fault-dip constants
    """
    C***** INPUT
    C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)
    C*****   DIP   : DIP-ANGLE (DEGREE)
    C### CAUTION ### IF COS(DIP) IS SUFFICIENTLY SMALL, IT IS SET TO ZERO
    """
    global ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD, S2D, C2D

    # Constants
    F0, F1, F2, PI2 = 0.0, 1.0, 2.0, 6.283185307179586
    EPS = 1.0e-6

    # Calculate constants
    ALP1 = (F1 - ALPHA) / F2
    ALP2 = ALPHA / F2
    ALP3 = (F1 - ALPHA) / ALPHA
    ALP4 = F1 - ALPHA
    ALP5 = ALPHA

    # Convert DIP angle to radians
    P18 = PI2 / 360.0
    SD = np.sin(DIP * P18)
    CD = np.cos(DIP * P18)

    # Handle the case where cos(DIP) is close to zero
    if np.abs(CD) < EPS:
        CD = F0
        if SD > F0:
            SD = F1
        elif SD < F0:
            SD = -F1

    SDSD = SD * SD
    CDCD = CD * CD
    SDCD = SD * CD
    S2D = F2 * SDCD
    C2D = CDCD - SDSD

    # Sets global variables ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD, S2D, C2D
    return


def DCCON1(X, Y, D):
    # Function to calculate station geometry constants for a point source
    global P, Q, S, T, XY, X2, Y2, D2, R, R2, R3, R5, R7, A3, A5, B3, C3, QR, QRX, UY, VY, WY, UZ, VZ, WZ

    # Constants
    F0, F1, F3, F5, EPS = 0.0, 1.0, 3.0, 5.0, 1.0e-6

    # Handle small values
    if np.abs(X) < EPS:
        X = F0
    if np.abs(Y) < EPS:
        Y = F0
    if np.abs(D) < EPS:
        D = F0

    P = Y * CD + D * SD
    Q = Y * SD - D * CD
    S = P * SD + Q * CD
    T = P * CD - Q * SD
    XY = X * Y
    X2 = X * X
    Y2 = Y * Y
    D2 = D * D
    R2 = X2 + Y2 + D2
    R = np.sqrt(R2)

    if R == F0:
        return

    R3 = R * R2
    R5 = R3 * R2
    R7 = R5 * R2

    A3 = F1 - F3 * X2 / R2
    A5 = F1 - F5 * X2 / R2
    B3 = F1 - F3 * Y2 / R2
    C3 = F1 - F3 * D2 / R2

    QR = F3 * Q / R5
    QRX = F5 * QR * X / R2

    UY = SD - F5 * Y * Q / R2
    UZ = CD + F5 * D * Q / R2
    VY = S - F5 * Y * P * Q / R2
    VZ = T + F5 * D * P * Q / R2
    WY = UY + SD
    WZ = UZ + CD

    # Sets global variables P, Q, S, T, XY, X2, Y2, D2, R, R2, R3, R5, R7, A3, A5, B3, C3, QR, QRX,
    # UY, VY, WY, UZ, VZ, WZ = [0.0] * 25
    return


def DC3D0(ALPHA, X, Y, Z, DEPTH, DIP, POT1, POT2, POT3, POT4):
    """
    C***** INPUT
    C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)
    C*****   X,Y,Z : COORDINATE OF OBSERVING POINT
    C*****   DEPTH : SOURCE DEPTH
    C*****   DIP   : DIP-ANGLE (DEGREE)
    C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY
    C*****       POTENCY=(  MOMENT OF DOUBLE-COUPLE  )/MYU     FOR POT1,2
    C*****       POTENCY=(INTENSITY OF ISOTROPIC PART)/LAMBDA  FOR POT3
    C*****       POTENCY=(INTENSITY OF LINEAR DIPOLE )/MYU     FOR POT4
    C
    C***** OUTPUT
    C*****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF POTENCY) /
    C*****               :                     (UNIT OF X,Y,Z,DEPTH)**2  )
    C*****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT= UNIT OF POTENCY) /
    C*****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH)**3  )
    C*****   UXZ,UYZ,UZZ : Z-DERIVATIVE
    C*****   IRET        : RETURN CODE
    C*****               :   =0....NORMAL
    C*****               :   =1....SINGULAR
    C*****               :   =2....POSITIVE Z WAS GIVEN
    """
    # Function to calculate displacement and strain due to a buried point source

    # Constants
    F0 = 0.0
    R = 0.0

    # Arrays
    U = np.zeros(12)

    # Return code
    IRET = 0

    # Check if Z is positive
    if Z > 0.:
        IRET = 2
        return IRET, np.zeros((1, 3)), np.zeros((3, 3))
    initialize_globals()

    # Real-source contribution
    XX, YY, ZZ, DD = X, Y, Z, DEPTH + Z

    # Call functions for real-source contribution
    DCCON0(ALPHA, DIP)
    DCCON1(XX, YY, DD)

    # Check if R is zero
    if R == F0:
        IRET = 1
        return IRET, np.zeros((1, 3)), np.zeros((3, 3))

    # Potencies
    PP1, PP2, PP3, PP4 = POT1, POT2, POT3, POT4

    # Calculate displacement for real-source contribution
    DUA = UA0(XX, YY, DD, PP1, PP2, PP3, PP4)
    for i in range(12):
        if i < 10:
            U[i] = U[i] - DUA[i]
        else:
            U[i] = U[i] + DUA[i]

    # Image-source contribution
    DD = DEPTH - Z
    DCCON1(XX, YY, DD)
    DUA = UA0(XX, YY, DD, PP1, PP2, PP3, PP4)
    DUB = UB0(XX, YY, DD, ZZ, PP1, PP2, PP3, PP4)
    DUC = UC0(XX, YY, DD, ZZ, PP1, PP2, PP3, PP4)

    # Calculate displacement for image-source contribution
    for i in range(12):
        DU = DUA[i] + DUB[i] + ZZ * DUC[i]
        if i >= 10:
            DU = DU + DUC[i - 9]
        U[i] = U[i] + DU

    # Assign values to variables
    u = np.array([U[0], U[1], U[2]])
    grad_u = np.array([U[3], U[4], U[5]], [U[6], U[7], U[8]], [U[9], U[10], U[11]])  # strain tensor components

    return IRET, u, grad_u
