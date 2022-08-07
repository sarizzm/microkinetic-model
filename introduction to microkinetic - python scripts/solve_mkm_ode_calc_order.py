import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt

#
# Python script to solve bimolecular reactions
# Section 4.2 - Bimolecular reaction
#


def main():
    """
    Calculate the reaction order in CO

    """
    trange = np.linspace(400, 1000, 25)
    # print(trange)
    T = 800 # temperature in K
    p,r = calc_order(T)
    m,b = np.polyfit(p, r, 1)
    plt.plot(p, r, 'o', label='Data points')
    plt.plot(p, m*p+b, '--', label='Linear fit')
    plt.legend()
    plt.xlabel('Pressure in Pa')
    plt.ylabel('log(rate)')
    plt.title('Reaction order CO = %f' % (m*p[2]))
    plt.show()

    ro_t =[]
    for tem in trange:
        p, r = calc_order(tem)
        m, b = np.polyfit(p, r, 1)
        ro_t.append(m * p[2])
    plt.plot(trange, ro_t, '-o', label='CO')
    # plt.plot(p, m * p + b, '--', label='Linear fit')
    plt.legend()
    plt.xlabel('Temperate in K')
    plt.ylabel('Reaction order')
    # plt.title('Reaction order CO = %f' % (m * p[2]))
    plt.show()


def calc_order(T):
    """
    18 Calculate reaction order at temperature T
    19 """
    pt = 20
    pa = 2.0/3.0 * pt * 1e5 # pressure of CO in Pa
    pb = 1.0/3.0 * pt * 1e5 # pressure of O2 in Pa
    pc = 0
    mc = 80 * 1.66054e-27

    # set series of factors to expend pressure in
    diffs = [0.95, 0.98, 1.0, 1.02, 1.05]
    rates = []
    for diff in diffs:
        x, y = solve_odes(T, pa * diff, pb, pc)
        r_co2 = calc_kdes(T, 1e-20, mc, 1, 0.561, 10e3) * y[-1,2]
        rates.append(r_co2)
    return np.multiply(diffs, pa), np.log(rates)




# def calc_range(T1, T2):
#     """
#     Calculate reaction rate over a range of temperatures
#     """
#     trange = np.linspace(T1, T2, 40)
#     rates = []
#     for T in trange:
#         x, y = solve_odes(T)
#         r_co2 = calc_kdes(T, 1e-20, 80 * 1.66054e-27, 1, 0.561, 10e3) * y[-1,2]
#         r_co =  calc_k_arr(T, 1e13, 120e3)* y[-1,0] * y[-1,1]
#         r_o2 = 0.5 * r_co
#         # print(y,'\n #################')
#
#         rates.append([r_co2,r_co,r_o2,y[-1,0],y[-1,1],y[-1,2],y[-1,3]])
#         # rates.append(r_co)
#         # rates.append(r_o2)
#
#
#     return trange, np.array(rates)
#


def solve_odes(T,pa, pb, pc):
    """
    Time-integrate chemo-kinetic system
    """
    # initial conditions
    y0 = [0, 0, 0, 1]
    t0 = 0
    t1 = 1   # total integration time
    pt = 20     # total pressure in bar
    # pa = 2.0/3.0 * pt * 1e5      # pressure of CO in Pa
    # pb = 1.0/3.0 * pt * 1e5      # pressure of O2 in Pa
    # pc = 0

    # construct ODE solver
    r = ode(dydt).set_integrator('vode', method='bdf',
           atol=1e-8, rtol=1e-8, nsteps=1000, with_jacobian=True)
    r.set_initial_value(y0, t0).set_f_params([T, pa, pb, pc])

    # integrate on a logaritmic scale
    xx = np.linspace(-12.0, np.log10(t1), int((np.log10(t1) + 12.0) * 10))
    # print(int((np.log10(t1) + 12.0) * 10))
    yy = []
    tt = []
    for x in xx:
        tnew = 10.0**x
        tt.append(tnew)
        yy.append(r.integrate(tnew))
    # print('yy is ', yy)

    return tt, np.matrix(yy)

def dydt(t, y, params):
    """
    Set of ordinary differential equations
    """
    T =  params[0]
    pa = params[1]
    pb = params[2]
    pc = params[3]

    dydt = np.zeros(4)

    ma = 28 * 1.66054e-27
    mb = 32 * 1.66054e-27
    mc = 80 * 1.66054e-27

    # calculate all reaction rate constants
    k_ads_1 = calc_kads(T, pa, 1e-20, ma)
    k_des_1 = calc_kdes(T, 1e-20, ma, 1, 2.8, 80e3)
    k_ads_2 = calc_kads(T, pb, 1e-20, mb)
    k_des_2 = calc_kdes(T, 1e-20, mb, 2, 2.08, 40e3)
    kf = calc_k_arr(T, 1e13, 120e3)
    kb = calc_k_arr(T, 1e13,  80e3)
    k_ads_4 = calc_kads(T, pc, 1e-20, mc)
    k_des_4 = calc_kdes(T, 1e-20, mc, 1, 0.561, 10e3)

    # collect similar terms in new variables
    r1f = k_ads_1 * y[3]
    r1b = k_des_1 * y[0]
    r2f = k_ads_2 * y[3]
    r2b = k_des_2 * y[1]**2
    r3f = kf * y[0] * y[1]
    r3b = kb * y[2] * y[3]
    r4f = k_ads_4 * y[3]
    r4b = k_des_4 * y[2]

    dydt[0] = r1f - r1b - r3f + r3b
    dydt[1] = 2.0 * r2f - 2.0 * r2b - r3f + r3b
    dydt[2] = r3f - r3b + r4f - r4b
    dydt[3] = -r1f + r1b - 2.0 * r2f + 2.0 * r2b + r3f - r3b - r4f + r4b

    return dydt

def calc_k_arr(T, nu, Eact):
    """
    Calculate reaction rate constant for a surface reaction

    T       - Temperature in K
    nu      - Pre-exponential factor in s^-1
    Eact    - Activation energy in J/mol
    """
    R = 8.3144598 # gas constant
    return nu * np.exp(-Eact / (R * T))

def calc_kads(T, P, A, m):
    """
    Reaction rate constant for adsorption

    T           - Temperature in K
    P           - Pressure in Pa
    A           - Surface area in m^2
    m           - Mass of reactant in kg
    """
    kb = 1.38064852E-23 # boltzmann constant
    return P*A / np.sqrt(2 * np.pi * m * kb * T)

def calc_kdes(T, A, m, sigma, theta_rot, Edes):
    """
    Reaction rate constant for desorption

    T           - Temperature in K
    A           - Surface area in m^2
    m           - Mass of reactant in kg
    sigma       - Symmetry number
    theta_rot   - Rotational temperature in K
    Edes        - Desorption energy in J/mol
    """
    kb = 1.38064852e-23 # boltzmann constant
    h = 6.62607004e-34  # planck constant
    R = 8.3144598       # gas constant
    return kb * T**3 / h**3 * A * (2 * np.pi * m * kb) / \
        (sigma * theta_rot) * np.exp(-Edes / (R*T))

if __name__ == '__main__':
    main()