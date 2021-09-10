---
layout: post
title: "Coursework: Features of a 2D Ising Model simulated with Markov chain Monte Carlo methods"
mathjax: true
published: true
---

## Context

This post documents a piece of coursework I completed in April 2021 for the final year of my physics batchelors degree. The brief was to use Python (or C++) to build an outlined simulation, or analyse some given data, extract meaningful resulsts, then write a 3000 word paper describing and evaluating the process.

There was not total freedom to choose the subject of the work, but I had read about the Ising model previously and was eager to learn more so this suited me fine.

The paper itself is attached as a pdf below, and I have reproduced the abstract for this webpage. Meanwhile, all the code is available on my [GitHub](https://github.com/JakeSkelton/ising-model-MCMC), and there is an example listing from the main numerical module below.

## Abstract

In this paper, we present observations of a two-dimensional Ising model, simulated with the Metropolis algorithm, a Markov chain Monte Carlo method. Square lattices are used, with periodic boundary conditions and nearest-neighbour spin-spin interactions. Considerations of Markov chain convergence and time correlations are taken into account, and the variation of magnetisation and heat capacity with temperature is observed, for a range of lattice sizes and applied magnetic fields. Hysteretic phenomena in the Ising model under an applied field are also investigated. In addition, we make an estimate of the critical temperature of an infinite 2D lattice, $$2.263 \pm 0.013 \; J/k_B$$. Onsager's analytical result of $$2/\ln(1 + \sqrt{2}) \; J/k_B$$ differs from our estimate by only 0.2%, and lies within error.

### [Download PDF](/assets/ising_model.pdf)

## Metropolis-Hastings algorithm subroutine

```cython
# -*- coding: utf-8 -*-
import numpy as np
cimport cython
from libc.math cimport exp

@cython.boundscheck(False)
@cython.nonecheck(False)
def metrohaste_vect(int numsteps, int[:,:] s, double[:] H, double[:] T, RNG):
    """
    A fast, Cython-based method to carry out one timestep of the Metropolis-
    -Hastings algorithm, a Markov chain Monte Carlo method. The lattice
    is stepped through in a left-to-right scanning fashion. This version
    outputs statistics of the lattice over time, like metrohaste_stats, but
    takes temperature and magnetic field arguments as vectors of length
    'numsteps', in order to cycle the magnetic field or anneal the lattice
    according to a predefined trajectory.

    Parameters
    ----------
    numsteps : int
        The number of Metropolis-Hastings timesteps to perform.
    s : square numpy array of ints
        The spin lattice, must be square and of type int.
    H : float
        The applied magnetic field, with dimensions (energy)/(spin).
    T : float
        The heat bath temperature, with dimensions of (energy).
    RNG : numpy.random generator object
        Seeded random number generator. Passing this to the function allows
        results to be repeated starting from a mother seed.

    Returns
    -------
    s_new : numpy array
        The updated spin lattice, ater 'numsteps' timesteps.
    sbars : numpy array    
        A 1D array of the lattice's mean spin at each timestep.
    energies : 
        A 1D array of the lattice's total energy at each timestep.
    """
    cdef Py_ssize_t N_ind = s.shape[0]
    cdef int N = s.shape[0]
    cdef double deltaE, E_t = 0, T_t = T[0], H_t = T[0]
    cdef int count = 0, sqcount = 0
    cdef double[:,:] p = np.zeros((N, N))
    cdef double[:] sbars = np.zeros(numsteps, dtype=float)
    cdef double[:] energy = np.zeros(numsteps, dtype=float)
    cdef Py_ssize_t t, i, j
    
    for i in range(N_ind):
        for j in range(N_ind):
            E_t -= (0.25*(s[(i-1)%N, j] + s[(i+1)%N, j] +
                    s[i, (j-1)%N] + s[i, (j+1)%N]) + 2*H_t)*s[i,j]
    energy[0] = E_t
    
    for t in range(numsteps):
        T_t = T[t]
        H_t = H[t]
        count = 0
        sqcount = 0
        p = RNG.uniform(0, 1, size=(N,N))
        for i in range(N_ind):
            for j in range(N_ind):
                deltaE = (s[(i-1)%N, j] + s[(i+1)%N, j] +
                          s[i, (j-1)%N] + s[i, (j+1)%N] + 2*H_t)*s[i,j]
                if (exp(-deltaE/T_t) > p[i,j]):
                    s[i,j] *= -1
                    E_t += deltaE
                count += s[i,j]
                sqcount += s[i,j]*s[i,j]
        sbars[t] = <double>count/<double>(N*N)
        energy[t] = E_t
        
    return (np.asarray(s), np.asarray(sbars), np.asarray(energy))
```

## A word on copyright

With the obvious exception of quoted material in the text, the paper, all it's figures, and the code are entirely my own work. You are welcome to re-use the figures (and of course the programs can generate more) so long as the source is attributed, and as GitHub will tell you, the code is under the GPL. 