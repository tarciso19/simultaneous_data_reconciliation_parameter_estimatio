# package importing
from motor import Rec
from casadi import vertcat

# model definition
def model(xr, param):

    T = vertcat(xr[6], xr[7], xr[8], xr[9], xr[10], xr[11], xr[25], xr[26], xr[27], xr[28], xr[29], xr[30],
                xr[44], xr[45], xr[46], xr[47], xr[48], xr[49], xr[63], xr[64], xr[65], xr[66], xr[67], xr[68],
                xr[82], xr[83], xr[84], xr[85], xr[86], xr[87])

    A, B, C, D , E = param[0], param[1], param[2], param[3], param[4]

    return A*T**2.5 + B*T**2 + C*T**1.5 + D*T + E

# process constraints definition
def constraints(xr, param):

    Fr1A, Fr2A, Fr3A, Fr4A, Fr5A, Fr6A = xr[0],xr[1],xr[2],xr[3],xr[4],xr[5]
    Tr1A, Tr2A, Tr3A, Tr4A, Tr5A, Tr6A = xr[6],xr[7],xr[8],xr[9],xr[10],xr[11]
    cr1A, cr2A, cr3A, cr4A, cr5A, cr6A = xr[12],xr[13],xr[14],xr[15],xr[16],xr[17]
    QrA = xr[18]
    Fr1B, Fr2B, Fr3B, Fr4B, Fr5B, Fr6B = xr[19], xr[20], xr[21], xr[22], xr[23], xr[24]
    Tr1B, Tr2B, Tr3B, Tr4B, Tr5B, Tr6B = xr[25], xr[26], xr[27], xr[28], xr[29], xr[30]
    cr1B, cr2B, cr3B, cr4B, cr5B, cr6B = xr[31], xr[32], xr[33], xr[34], xr[35], xr[36]
    QrB = xr[37]
    Fr1C, Fr2C, Fr3C, Fr4C, Fr5C, Fr6C = xr[38], xr[39], xr[40], xr[41], xr[42], xr[43]
    Tr1C, Tr2C, Tr3C, Tr4C, Tr5C, Tr6C = xr[44], xr[45], xr[46], xr[47], xr[48], xr[49]
    cr1C, cr2C, cr3C, cr4C, cr5C, cr6C = xr[50], xr[51], xr[52], xr[53], xr[54], xr[55]
    QrC = xr[56]
    Fr1D, Fr2D, Fr3D, Fr4D, Fr5D, Fr6D = xr[57], xr[58], xr[59], xr[60], xr[61], xr[62]
    Tr1D, Tr2D, Tr3D, Tr4D, Tr5D, Tr6D = xr[63], xr[64], xr[65], xr[66], xr[67], xr[68]
    cr1D, cr2D, cr3D, cr4D, cr5D, cr6D = xr[69], xr[70], xr[71], xr[72], xr[73], xr[74]
    QrD = xr[75]
    Fr1E, Fr2E, Fr3E, Fr4E, Fr5E, Fr6E = xr[76], xr[77], xr[78], xr[79], xr[80], xr[81]
    Tr1E, Tr2E, Tr3E, Tr4E, Tr5E, Tr6E = xr[82], xr[83], xr[84], xr[85], xr[86], xr[87]
    cr1E, cr2E, cr3E, cr4E, cr5E, cr6E = xr[88], xr[89], xr[90], xr[91], xr[92], xr[93]
    QrE = xr[94]

    A, B, C, D, E = param[0], param[1], param[2], param[3], param[4]

    return vertcat(Fr1A - Fr2A - Fr3A,
                   Fr2A - Fr4A,
                   Fr3A - Fr5A,
                   Fr4A + Fr5A - Fr6A,
                   Tr1A - Tr2A,
                   cr1A*Fr1A*(Tr1A-25) - cr2A*Fr2A*(Tr2A-25) - cr3A*Fr3A*(Tr3A-25),
                   cr2A*Fr2A*(Tr2A-25) - cr4A*Fr4A*(Tr4A-25) + QrA*1000,
                   cr3A*Fr3A*(Tr3A-25) - cr5A*Fr5A*(Tr5A-25),
                   cr4A*Fr4A*(Tr4A-25) + cr5A*Fr5A*(Tr5A-25) - cr6A*Fr6A*(Tr6A-25),
                   A*Tr1A**2.5 + B*Tr1A**2 + C*Tr1A**1.5 + D*Tr1A + E - cr1A,
                   A*Tr2A**2.5 + B*Tr2A**2 + C*Tr2A**1.5 + D*Tr2A + E - cr2A,
                   A*Tr3A**2.5 + B*Tr3A**2 + C*Tr3A**1.5 + D*Tr3A + E - cr3A,
                   A*Tr4A**2.5 + B*Tr4A**2 + C*Tr4A**1.5 + D*Tr4A + E - cr4A,
                   A*Tr5A**2.5 + B*Tr5A**2 + C*Tr5A**1.5 + D*Tr5A + E - cr5A,
                   A*Tr6A**2.5 + B*Tr6A**2 + C*Tr6A**1.5 + D*Tr6A + E - cr6A,
                   Fr1B - Fr2B - Fr3B,
                   Fr2B - Fr4B,
                   Fr3B - Fr5B,
                   Fr4B + Fr5B - Fr6B,
                   Tr1B - Tr2B,
                   cr1B * Fr1B * (Tr1B - 25) - cr2B * Fr2B * (Tr2B - 25) - cr3B * Fr3B * (Tr3B - 25),
                   cr2B * Fr2B * (Tr2B - 25) - cr4B * Fr4B * (Tr4B - 25) + QrB*1000,
                   cr3B * Fr3B * (Tr3B - 25) - cr5B * Fr5B * (Tr5B - 25),
                   cr4B * Fr4B * (Tr4B - 25) + cr5B * Fr5B * (Tr5B - 25) - cr6B * Fr6B * (Tr6B - 25),
                   A * Tr1B ** 2.5 + B * Tr1B ** 2 + C * Tr1B ** 1.5 + D * Tr1B + E - cr1B,
                   A * Tr2B ** 2.5 + B * Tr2B ** 2 + C * Tr2B ** 1.5 + D * Tr2B + E - cr2B,
                   A * Tr3B ** 2.5 + B * Tr3B ** 2 + C * Tr3B ** 1.5 + D * Tr3B + E - cr3B,
                   A * Tr4B ** 2.5 + B * Tr4B ** 2 + C * Tr4B ** 1.5 + D * Tr4B + E - cr4B,
                   A * Tr5B ** 2.5 + B * Tr5B ** 2 + C * Tr5B ** 1.5 + D * Tr5B + E - cr5B,
                   A * Tr6B ** 2.5 + B * Tr6B ** 2 + C * Tr6B ** 1.5 + D * Tr6B + E - cr6B,
                   Fr1C - Fr2C - Fr3C,
                   Fr2C - Fr4C,
                   Fr3C - Fr5C,
                   Fr4C + Fr5C - Fr6C,
                   Tr1C - Tr2C,
                   cr1C*Fr1C*(Tr1C-25) - cr2C*Fr2C*(Tr2C-25) - cr3C*Fr3C*(Tr3C-25),
                   cr2C*Fr2C*(Tr2C-25) - cr4C*Fr4C*(Tr4C-25) + QrC*1000,
                   cr3C*Fr3C*(Tr3C-25) - cr5C*Fr5C*(Tr5C-25),
                   cr4C*Fr4C*(Tr4C-25) + cr5C*Fr5C*(Tr5C-25) - cr6C*Fr6C*(Tr6C-25),
                   A*Tr1C**2.5 + B*Tr1C**2 + C*Tr1C**1.5 + D*Tr1C + E - cr1C,
                   A*Tr2C**2.5 + B*Tr2C**2 + C*Tr2C**1.5 + D*Tr2C + E - cr2C,
                   A*Tr3C**2.5 + B*Tr3C**2 + C*Tr3C**1.5 + D*Tr3C + E - cr3C,
                   A*Tr4C**2.5 + B*Tr4C**2 + C*Tr4C**1.5 + D*Tr4C + E - cr4C,
                   A*Tr5C**2.5 + B*Tr5C**2 + C*Tr5C**1.5 + D*Tr5C + E - cr5C,
                   A*Tr6C**2.5 + B*Tr6C**2 + C*Tr6C**1.5 + D*Tr6C + E - cr6C,
                   Fr1D - Fr2D - Fr3D,
                   Fr2D - Fr4D,
                   Fr3D - Fr5D,
                   Fr4D + Fr5D - Fr6D,
                   Tr1D - Tr2D,
                   cr1D * Fr1D * (Tr1D - 25) - cr2D * Fr2D * (Tr2D - 25) - cr3D * Fr3D * (Tr3D - 25),
                   cr2D * Fr2D * (Tr2D - 25) - cr4D * Fr4D * (Tr4D - 25) + QrD*1000,
                   cr3D * Fr3D * (Tr3D - 25) - cr5D * Fr5D * (Tr5D - 25),
                   cr4D * Fr4D * (Tr4D - 25) + cr5D * Fr5D * (Tr5D - 25) - cr6D * Fr6D * (Tr6D - 25),
                   A * Tr1D ** 2.5 + B * Tr1D ** 2 + C * Tr1D ** 1.5 + D * Tr1D + E - cr1D,
                   A * Tr2D ** 2.5 + B * Tr2D ** 2 + C * Tr2D ** 1.5 + D * Tr2D + E - cr2D,
                   A * Tr3D ** 2.5 + B * Tr3D ** 2 + C * Tr3D ** 1.5 + D * Tr3D + E - cr3D,
                   A * Tr4D ** 2.5 + B * Tr4D ** 2 + C * Tr4D ** 1.5 + D * Tr4D + E - cr4D,
                   A * Tr5D ** 2.5 + B * Tr5D ** 2 + C * Tr5D ** 1.5 + D * Tr5D + E - cr5D,
                   A * Tr6D ** 2.5 + B * Tr6D ** 2 + C * Tr6D ** 1.5 + D * Tr6D + E - cr6D,
                   Fr1E - Fr2E - Fr3E,
                   Fr2E - Fr4E,
                   Fr3E - Fr5E,
                   Fr4E + Fr5E - Fr6E,
                   Tr1E - Tr2E,
                   cr1E * Fr1E * (Tr1E - 25) - cr2E * Fr2E * (Tr2E - 25) - cr3E * Fr3E * (Tr3E - 25),
                   cr2E * Fr2E * (Tr2E - 25) - cr4E * Fr4E * (Tr4E - 25) + QrE*1000,
                   cr3E * Fr3E * (Tr3E - 25) - cr5E * Fr5E * (Tr5E - 25),
                   cr4E * Fr4E * (Tr4E - 25) + cr5E * Fr5E * (Tr5E - 25) - cr6E * Fr6E * (Tr6E - 25),
                   A * Tr1E ** 2.5 + B * Tr1E ** 2 + C * Tr1E ** 1.5 + D * Tr1E + E - cr1E,
                   A * Tr2E ** 2.5 + B * Tr2E ** 2 + C * Tr2E ** 1.5 + D * Tr2E + E - cr2E,
                   A * Tr3E ** 2.5 + B * Tr3E ** 2 + C * Tr3E ** 1.5 + D * Tr3E + E - cr3E,
                   A * Tr4E ** 2.5 + B * Tr4E ** 2 + C * Tr4E ** 1.5 + D * Tr4E + E - cr4E,
                   A * Tr5E ** 2.5 + B * Tr5E ** 2 + C * Tr5E ** 1.5 + D * Tr5E + E - cr5E,
                   A * Tr6E ** 2.5 + B * Tr6E ** 2 + C * Tr6E ** 1.5 + D * Tr6E + E - cr6E)


# initial estimative
x0 = [101.9,64.4,34.6,64.2,36.4,99,45.1,47.3,43,78.1,46.4,66.2,4.19,4.19,4.19,4.19,4.19,4.19,8,
      101.9,64.4,34.6,64.2,36.4,99,18.6,15.6,23.6,71.3,23.7,54.6,4.18,4.18,4.18,4.18,4.18,4.18,15,
      101.9,64.4,34.6,64.2,36.4,99,45.1,43.3,50,80,49,69,4.18,4.18,4.18,4.18,4.18,4.18,10,
      101.9,64.4,34.6,64.2,36.4,99,45.1,45,46.9,59.9,45.8,55.3,4.18,4.18,4.18,4.18,4.18,4.18,4,
      101.9,64.4,34.6,64.2,36.4,99,30.1,27.3,35,72.1,34.4,59.2,4.19,4.19,4.19,4.19,4.19,4.19,12,
      0,0,0,0,4]


# initializing the calculation engine
SDRPE = Rec(constraints=constraints, model=model)

# data definition
SDRPE.dataRead(file=r'C:\Users\Tarciso\Desktop\data.xlsx', sheet='DR', dataType='DR') # please Change the "file" argument according to
                                                                                   # the path where the "data.xlsx" archive is on your computer.
SDRPE.dataRead(file=r'C:\Users\Tarciso\Desktop\data.xlsx', sheet='PE', dataType='PE') # please Change the "file" argument according to
                                                                                   # the path where the "data.xlsx" archive is on your computer.
# SDRPE resolution
SDRPE.simultaneousDataReconciliationParameterEstimation(initial_estimative=x0, lbg=[0]*15*5, ubg=[0]*15*5, algoritmo='ipopt'
                                                        ,method='decoupled',first='PE')

# results report
SDRPE.report('Decoupled Method - First_Data Reconciliation')
# graph results
SDRPE.charts()