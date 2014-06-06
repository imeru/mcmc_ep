

# Import relevant modules
import numpy as np
from scipy.stats import triang


def tri_logdensity(x, min, max, mode):
    loc = min
    scale = max - min
    c = (mode - min)/scale
    log_density = triang.logpdf(x = x, c = c, loc = loc, scale = scale)
    return log_density


# Priors 
def prior(param):
    ROOF = param[0]
    WALL = param[1]
    EPD = param[2]
    LPD = param[3]
    HSP = param[4]
    INF = param[5]
    Boiler = param[6]
    COP = param[7]

    ROOF_prior = tri_logdensity(x = ROOF, min = 0.01, max = 0.25, mode = 0.09667)
    WALL_prior = tri_logdensity(x = WALL, min = 0.01, max = 0.25, mode = 0.055)
    EPD_prior = tri_logdensity(x = EPD, min = 1, max = 70, mode = 22.8)
    LPD_prior = tri_logdensity(x = LPD, min = 1, max = 50, mode = 14.5)
    HSP_prior = tri_logdensity(x = HSP, min = 17, max = 28, mode = 21)
    INF_prior = tri_logdensity(x = INF, min = 0.1, max = 2, mode = 0.675)
    Boiler_prior = tri_logdensity(x = Boiler, min = 0.5, max = 0.99, mode = 0.72)
    COP_prior = tri_logdensity(x = COP, min = 2, max = 4, mode = 2.65)
    return ROOF_prior+WALL_prior+EPD_prior+LPD_prior+HSP_prior+INF_prior+Boiler_prior+COP_prior



# likelihoods

def likelihooss(param):
ROOF = param[0]
