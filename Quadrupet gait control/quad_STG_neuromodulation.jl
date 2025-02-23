#=
This file generates an initial set of STG neurons using the STG model as well as
neuromodulating this set to tune its firing pattern
=#

# Include DICs computation
include("quad_STG_DIC.jl") # Include the computation of DICs and compensation algorithms

## Functions to generate the initial set of STG neurons
# This function create 'ncells' STG neurons with a maximal variability in terms of
# maximal conductances while maintaining a stable firing pattern (here spiking)
function degeneracy_fixDICs_neuromod(ncells, gsth, guth, Vth, tmKCa)
  # Initialiizing some variables
  g_all_init = zeros(ncells, 8)
  ICs_th_init = zeros(ncells, 8)

  # Fitted linear relation
  gfth = - gsth - 2.2

  # Loop over all STG neurons
  for n = 1 : ncells
    # Generating random gleak and gCaT, gCaS, gKd, gKCa (proportional to gleak)
    gleak = 0.007 + 0.007 * rand(1, 1)[1]
    gCaT = gleak * (2. + 5. * rand(1, 1)[1]) / 0.01
    gCaS = gleak * (6. + 16. * rand(1, 1)[1]) / 0.01
    gKd = gleak * (140. + 40. * rand(1, 1)[1]) / 0.01
    gKCa = gleak * (70. + 70. * rand(1, 1)[1]) / 0.01

    # Computing gNa, gA and gH using the compensation algorithm to observe stable spiking
    (gNa, gA, gH) = DICs_gmax_initNaAH(gCaT, gCaS, gKd, gKCa, gleak, gfth, gsth, guth, Vth, tmKCa[n])

    # Computing DICs & co of the current STG neuron
    (ith, iosc, gf, gs, gu, gin, Istatic) = DICs(V, 0., gNa, gCaT, gCaS, gA, gKCa, gKd, gH, gleak, tmKCa[n])

    # Saving DICs & co at threshold and up state voltage and all maximal conductances
    ICs_th_init[n, :] = [V[ith], V[iosc], gf[ith], gs[ith], gs[iosc], gu[ith], gin[ith], Istatic[ith]]
    g_all_init[n, :] = [gNa, gCaT, gCaS, gA, gKCa, gKd, gH, gleak]
  end

  return g_all_init, ICs_th_init
end

## This function neuromodulate a set of STG neurons to switch from a certain firing pattern
# to another in a robust way using the compensation algorithm with gCaS and gA
function neuromodCaSA(ncells, g_all, ICs, gsth, guth, tmKCa)
  # Initialiizing some variables
  g_all_neuromod = zeros(ncells, 8)
  ICs_th_neuromod = zeros(ncells, 8)

  # Loop over all STG neurons
  for n = 1 : ncells
    # Extracting fixed maximal conductances as well as threshold voltage
    gNa = g_all[n, 1]
    gCaT = g_all[n, 2]
    gKCa = g_all[n, 5]
    gKd = g_all[n, 6]
    gH = g_all[n, 7]
    gleak = g_all[n, 8]
    Vth = ICs[n, 1]

    # Computing new gCaS and gA using the compensation algorithm to neuromodulate firing pattern
    (gCaS, gA) = DICs_gmax_neuromodCaSA(gNa, gCaT, gKd, gKCa, gH, gleak, gsth, guth, Vth, tmKCa[n])

    # Computing DICs & co of the current STG neuron
    (ith, iosc, gf, gs, gu, gin, Istatic) = DICs(V, 0., gNa, gCaT, gCaS, gA, gKCa, gKd, gH, gleak, tmKCa[n])

    # Saving DICs & co at threshold and up state voltage and all maximal conductances
    ICs_th_neuromod[n, :] = [V[ith], V[iosc], gf[ith], gs[ith], gs[iosc], gu[ith], gin[ith], Istatic[ith]]
    g_all_neuromod[n, :] = [gNa, gCaT, gCaS, gA, gKCa, gKd, gH, gleak]
  end

  return g_all_neuromod, ICs_th_neuromod
end
