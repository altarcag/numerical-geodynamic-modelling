import numpy as np

def Gibbs_MgO(P, T):
    R = 8.314
    Pr = 100000
    Tr = 298.15
    Hr = -601500.00
    Vr = 1.12228e-5
    phi = 30179500000
    c = [1.96612, 4.12756, 0.53690]
    dH = [2966.88, 5621.69, 27787.19]
    dV = [3.52971e-8, 3.52971e-8, 1.9849568e-6]

    F = 5/4 * (Pr + phi)**0.2 * ((P + phi)**0.8 - (Pr + phi)**0.8)
    Gm = Hr + Vr * F
    for i in range(3):
        e = np.exp(-(dH[i] + dV[i] * F) / (R*T))
        e0 = np.exp(-dH[i] / (R*Tr))
        Gm += c[i] * (R*T *np.log(1 - e) - dH[i]*e0 / (1 - e0))
    return Gm