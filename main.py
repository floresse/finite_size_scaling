from ctypes import *

simulation = CDLL('src/simulation.so')


def main():
    L = 64
    sigma = c_double(0.1)

    T_ini = c_double(10.7)
    T_fin = T_ini
    dT = c_double(0.1)

    h_ini = c_double(0.0)
    h_fin = h_ini
    dh = c_double(0.1)

    simulation.simulation(L, sigma, T_ini, T_fin, dT, h_ini, h_fin, dh)


if __name__ == '__main__':
    main()
