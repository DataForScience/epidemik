
from epidemik import EpiModel
import numpy as np

beta = 0.3
mu = 0.1

SIR0 = EpiModel()
SIR0.add_interaction('S', 'I', 'I', beta)
SIR0.add_spontaneous('I', 'R', mu)


SIR = EpiModel()
SIR.add_interaction('S', 'I', 'I', beta)
SIR.add_spontaneous('I', 'R', mu)

SIR2 = EpiModel()
SIR2.add_interaction('S', 'I', 'I', beta)
SIR2.add_spontaneous('I', 'R', mu)

print(SIR)

def make_seasonality(a_max=1.1, a_min=0.1, t_max=200):
    t = np.arange(1, 366)
    return a_min+1/2*(1+np.cos(2*np.pi/365*(t-t_max)))*(a_max-a_min)


season = make_seasonality()

def seasonality(t, population, pos):
	return season[t]

SIR0.integrate(365, S=10000, I=10)

SIR0.plot()


SIR.integrate(365, S=10000, I=10, seasonality=season)

SIR.plot()

SIR2.integrate(365, S=10000, I=10, seasonality=seasonality)

SIR2.plot()