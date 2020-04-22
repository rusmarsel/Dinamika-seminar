from __future__ import division
from numpy import linspace, array, append
import matplotlib.pyplot as plt

def func(t, M0=1, T=1):
    if t < (T / 2):
        return (2*t / T) * M0
    else:
        return (2*t / T) * M0 - 2*M0

if __name__ == "__main__":
    time_end = 10
    max_amp = 10
    t = linspace(0, time_end, 100)
    Mt = array([])
    for i in range(t.size):
        y = func(t[i], M0=max_amp, T=time_end)
        Mt = append(Mt, y)
    
    plt.figure()
    plt.plot(t, Mt)
    plt.show()
