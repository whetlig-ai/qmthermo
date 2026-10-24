# qmthermo
## Program description
qmthermo Is a Python 3 library, created to apply correction schemes accounting for the low frequency entropy "explosion".
The list of available models:
- HO (classic harmonic oscillator approximation)
- sRRHO (scaled model of ideal-gas rigid-rotor harmonic-oscillator approximation with scaling factor defined by the used method)
- msRRHO(t) (modified sRRHO approximation including the interpolation of the vibrational and rotational terms, parameter `t` is the frequency, when the contribution of both terms is equal). 
