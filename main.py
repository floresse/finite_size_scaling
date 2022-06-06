from ctypes import *

simulation = CDLL('src/simulation.so')


def main():
    L = 64
    sigma = c_double(0.1)

    T_ini = c_double(12.0)
    T_fin = c_double(28.0)
    dT = c_double(0.5)

    h_ini = c_double(0.0)
    h_fin = h_ini
    dh = c_double(0.1)

    simulation.simulation(L, sigma, T_ini, T_fin, dT, h_ini, h_fin, dh)


if __name__ == '__main__':
    main()
