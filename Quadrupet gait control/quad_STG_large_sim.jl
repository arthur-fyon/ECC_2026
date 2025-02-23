using Plots, LaTeXStrings, Random, Dierckx, DelimitedFiles
include("quad_STG_kinetics.jl") # Loading of STG kinetics of gating variables
include("quad_STG_models.jl") # Loading of STG model
include("quad_STG_utils.jl") # Loading of some utils functions
include("quad_STG_gs_derivatives.jl") # Loading of X_inf derivatives
include("quad_STG_DIC.jl") # Loading of the DIC and compensation algorithm
include("quad_STG_neuromodulation.jl"); # Loading of the neuromodulation cells functions

# Definition of simulation time (in ms)
const Tfinal = 50000
const dt = 0.0005
const tsim = 0 : dt : Tfinal
const dtplot = 0.2
const tt = 0 : dtplot : Tfinal
const dtratio = Int(dtplot/dt)
const Tdt = Int(Tfinal/dt)
const tt_index = 1 : dtratio : Tdt+1

# Definition of reversal potential values (in mV) and membrane capacitance
const VNa   = 50. # Sodium reversal potential
const VK    = -80. # Potassium reversal potential
const VCa   = 80. # Calcium reversal potential
const VH    = -20. # Reversal potential for the H-current (permeable to both sodium and potassium ions)
const Vleak = -50. # Reversal potential of leak channels
const Vsyn = -75. # Reversal potential of synaptic channels
const C     = 1. # Membrane capacitance

# Definition of voltage range for the DICs
const Vmin = -60
const Vmax = 0
const V    = range(Vmin, stop=Vmax, step=0.01)

# Fixing random seed
for s = 1 : 199
    display(s)
    Random.seed!(s)

    # Initial firing pattern
    guth = 4.
    Vth = -50.
    tmKCavec = 4*ones(8)
    (g_all_init, ICs_th_init) = degeneracy_fixDICs_neuromod(1, 5., guth, Vth, tmKCavec)
    # create a spiking set with max variability in gCaS and gA

    # Definition of parameters
    Iappvec = 0. * ones(8)
    gNavec = g_all_init[1] * ones(8)
    gCaTvec = g_all_init[2] * ones(8)
    gCaSvec = g_all_init[3] * ones(8)
    gAvec = g_all_init[4] * ones(8)
    gKCavec = g_all_init[5] * ones(8)
    gKdvec = g_all_init[6] * ones(8)
    gHvec = g_all_init[7] * ones(8)
    gleakvec = g_all_init[8] * ones(8)

    α = 5e-3 # Rate of transfer between intracellular and membrane
    β = 5e-3 # Rate of degradation of intracellular proteins
    Kp = 3e-4 # Proprtional gain
    Ki = 5e-6 # Integral gain
    Kt = β / Ki # Anti-windup gain
    gsth_sim_burst(t) = -8.
    gsth_sim_spike(t) = 5.
    gsth_spike2burst(t) = 5. - 13. * (t>Tfinal/2)
    gsth_burst2spike(t) = -8. + 13. * (t>Tfinal/2)
    gsthvec_trot2gallop = [gsth_sim_burst, gsth_sim_burst, gsth_sim_burst, gsth_sim_burst,
                           gsth_spike2burst, gsth_spike2burst, gsth_burst2spike, gsth_burst2spike]
    gsthvec_gallop2trot = [gsth_sim_burst, gsth_sim_burst, gsth_sim_burst, gsth_sim_burst,
                           gsth_burst2spike, gsth_burst2spike, gsth_spike2burst, gsth_spike2burst]
    guth_sim(t) = 4.
    guthvec = [guth_sim, guth_sim, guth_sim, guth_sim, guth_sim, guth_sim, guth_sim, guth_sim]
    u_maxCaS = 1e7 # Maximum value of actuator
    u_maxA = 1e7
    Vthvec = ICs_th_init[1] * ones(8)

    gsyn_HC = 0.8
    gsyn_drive = 0.8
    gsyn = [gsyn_HC; gsyn_HC; gsyn_HC; gsyn_HC; gsyn_HC; gsyn_HC; gsyn_HC; gsyn_HC;
            gsyn_drive; gsyn_drive; gsyn_drive; gsyn_drive; gsyn_drive; gsyn_drive; gsyn_drive; gsyn_drive]


    @time Vgallop2trot = simulateSTG_network(Iappvec, gNavec, gCaTvec, gCaSvec, gAvec, gKCavec, gKdvec,
                                              gHvec, gleakvec, gsyn, tmKCavec, α, β, Kp, Ki, Kt,
                                              gsthvec_gallop2trot, guthvec, Vthvec, u_maxCaS, u_maxA);

    spike_times_firstLH = extract_firstspike_times(Vgallop2trot[tt_index, 2], tt)
    spike_times_firstLF = extract_firstspike_times(Vgallop2trot[tt_index, 1], tt)
    spike_times_firstRH = extract_firstspike_times(Vgallop2trot[tt_index, 4], tt)
    spike_times_firstRF = extract_firstspike_times(Vgallop2trot[tt_index, 3], tt)

    phase1, t1, phase2, t2, phase3, t3, phase4, t4 = computing_phases(spike_times_firstLH,
        spike_times_firstLF, spike_times_firstRH, spike_times_firstRF)

    i1 = findall(t1 .> 15000 .&& t1 .< 45000)
    i2 = findall(t2 .> 15000 .&& t2 .< 45000)
    i3 = findall(t3 .> 15000 .&& t3 .< 45000)
    i4 = findall(t4 .> 15000 .&& t4 .< 45000)

    writedlm("./data/phase1_"*string(s)*".dat", phase1[i1])
    writedlm("./data/phase2_"*string(s)*".dat", phase2[i2])
    writedlm("./data/phase3_"*string(s)*".dat", phase3[i3])
    writedlm("./data/phase4_"*string(s)*".dat", phase4[i4])
    writedlm("./data/t1_"*string(s)*".dat", t1[i1])
    writedlm("./data/t2_"*string(s)*".dat", t2[i2])
    writedlm("./data/t3_"*string(s)*".dat", t3[i3])
    writedlm("./data/t4_"*string(s)*".dat", t4[i4])
    writedlm("./data/g_"*string(s)*".dat", g_all_init)
end
