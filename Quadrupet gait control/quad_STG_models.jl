#=
This file contains differential equations describing the STG model
=#

include("quad_STG_kinetics.jl") # Include STG model gating functions

# STG ODEs
function dV(V, mNa, hNa, mCaT, hCaT, mCaS, hCaS, mA, hA, mKCa, mKd, mH,
            Ca, Iapp, gNa, gCaT, gCaS, gA, gKCa, gKd, gH, gleak)
  (dt) * (1/C) * (- gNa*mNa^3*hNa*(V-VNa) - gCaT*mCaT^3*hCaT*(V-VCa) -
                    gCaS*mCaS^3*hCaS*(V-VCa) - gA*mA^3*hA*(V-VK) -
                    gKCa*mKCa^4*(V-VK) - gKd*mKd^4*(V-VK) - gH*mH*(V-VH) -
                    gleak*(V-Vleak) + Iapp)
end
dmNa(V, mNa) = (dt) * ((1/tau_mNa(V)) * (mNa_inf(V) - mNa))
dhNa(V, hNa) = (dt) * ((1/tau_hNa(V)) * (hNa_inf(V) - hNa))
dmCaT(V, mCaT) = (dt) * ((1/tau_mCaT(V)) * (mCaT_inf(V) - mCaT))
dhCaT(V, hCaT) = (dt) * ((1/tau_hCaT(V)) * (hCaT_inf(V) - hCaT))
dmCaS(V, mCaS) = (dt) * ((1/tau_mCaS(V)) * (mCaS_inf(V) - mCaS))
dhCaS(V, hCaS) = (dt) * ((1/tau_hCaS(V)) * (hCaS_inf(V) - hCaS))
dmA(V, mA) = (dt) * ((1/tau_mA(V)) * (mA_inf(V) - mA))
dhA(V, hA) = (dt) * ((1/tau_hA(V)) * (hA_inf(V) - hA))
dmKCa(V, Ca, mKCa, tmKCa) = (dt) * ((1/tau_mKCa(V, tmKCa)) * (mKCa_inf(V, Ca) - mKCa))
dmKd(V, mKd) = (dt) * ((1/tau_mKd(V)) * (mKd_inf(V) - mKd))
dmH(V, mH) = (dt) * ((1/tau_mH(V)) * (mH_inf(V) - mH))
dCa(V, mCaT, hCaT, mCaS, hCaS, Ca, gCaT, gCaS) = (dt) * (-0.94 * (gCaT*mCaT^3*hCaT*(V-VCa) + gCaS*mCaS^3*hCaS*(V-VCa)) - Ca + 0.05)/20
ds(V, s) = (dt) * ((1/tau_s(V)) * (s_inf(V) - s))
dgi(g, gi, α, β, u) = (dt) * (α * (g - gi) - β * gi + u)
dg(g, gi, α) = (dt) * (α * (gi - g))
dz(e, Kt, u, v) = (dt) * (e + Kt * (u - v))

# Control signal function
function uve(g, g_r, z, Kp, Ki, u_max)
    # Compute control signal
    e = g_r - g
    v = Kp * e + Ki * z

    # Anti-windup system
    if v > u_max
        u = u_max
    else
        u = v
    end

    return u, v, e
end

function simulateSTG_network(Iapp, gNa, gCaT, gCaSinit, gAinit, gKCa, gKd, gH, gleak, gsyn,
                             tmKCa, α, β, Kp, Ki, Kt, gsth, guth, Vth, u_maxCaS, u_maxA)

    # Extracting synaptic conductances
    gsyn12 = gsyn[1]
    gsyn21 = gsyn[2]
    gsyn34 = gsyn[3]
    gsyn43 = gsyn[4]
    gsyn56 = gsyn[5]
    gsyn65 = gsyn[6]
    gsyn78 = gsyn[7]
    gsyn87 = gsyn[8]
    gsyn51 = gsyn[9]
    gsyn53 = gsyn[10]
    gsyn62 = gsyn[11]
    gsyn64 = gsyn[12]
    gsyn72 = gsyn[13]
    gsyn73 = gsyn[14]
    gsyn81 = gsyn[15]
    gsyn84 = gsyn[16]

    # Initial conditions
    gCaS = copy(gCaSinit)
    gCaSi = copy(gCaSinit)
    zCaS = zero(gCaS)
    gA = copy(gAinit)
    gAi = copy(gAinit)
    zA = zero(gA)

    Vprev = .-60. .* ones(8) .+ 10. .* (rand(8).-0.5)
    Vcur = copy(Vprev)
    mNa = mNa_inf.(Vprev)
    hNa = hNa_inf.(Vprev)
    mCaT = mCaT_inf.(Vprev)
    hCaT = hCaT_inf.(Vprev)
    mCaS = mCaS_inf.(Vprev)
    hCaS = hCaS_inf.(Vprev)
    mA = mA_inf.(Vprev)
    hA = hA_inf.(Vprev)
    Ca = .-0.94 .* (gCaT .* mCaT.^3 .* hCaT .* (Vprev.-VCa) .+
                    gCaS .* mCaS.^3 .* hCaS .* (Vprev.-VCa)) .+ 0.05
    mKCa = mKCa_inf.(Vprev, Ca)
    mKd = mKd_inf.(Vprev)
    mH = mH_inf.(Vprev)

    gs_th = zeros(Tdt+1, 8)
    gu_th = zero(gs_th)
    for j = 1 : 8
        gs_th[:, j] = gsth[j].(tsim)
        gu_th[:, j] = guth[j].(tsim)
    end

    gCaS_r = zeros(8)
    gA_r = zeros(8)
    for j = 1 : 8
        (gCaS_r[j], gA_r[j]) = DICs_gmax_neuromodCaSA(gNa[j], gCaT[j], gKd[j], gKCa[j], gH[j],
                                                gleak[j], gs_th[1, j], gu_th[1, j], Vth[j], tmKCa[j])
    end

    s12 = s_inf(Vprev[1])
    s21 = s_inf(Vprev[2])
    s34 = s_inf(Vprev[3])
    s43 = s_inf(Vprev[4])
    s56 = s_inf(Vprev[5])
    s65 = s_inf(Vprev[6])
    s78 = s_inf(Vprev[7])
    s87 = s_inf(Vprev[8])
    s51 = s_inf(Vprev[5])
    s53 = s_inf(Vprev[5])
    s62 = s_inf(Vprev[6])
    s64 = s_inf(Vprev[6])
    s72 = s_inf(Vprev[7])
    s73 = s_inf(Vprev[7])
    s81 = s_inf(Vprev[8])
    s84 = s_inf(Vprev[8])

    # Initialize saving variables
    V_sol = zeros(Tdt+1, 8)
    V_sol[1, :] = copy(Vcur)

    for z = 2 : Tdt+1
        # STG ODEs
        for j = 1 : 8
            # Liu model ODEs
            Vcur[j] += dV(Vprev[j], mNa[j], hNa[j], mCaT[j], hCaT[j], mCaS[j], hCaS[j], mA[j], hA[j], mKCa[j], mKd[j],
                          mH[j], Ca[j], Iapp[j], gNa[j], gCaT[j], gCaS[j], gA[j], gKCa[j], gKd[j], gH[j], gleak[j])
            mKCa[j] += dmKCa(Vprev[j], Ca[j], mKCa[j], tmKCa[j])
            Ca[j] += dCa(Vprev[j], mCaT[j], hCaT[j], mCaS[j], hCaS[j], Ca[j], gCaT[j], gCaS[j])
            mNa[j] += dmNa(Vprev[j], mNa[j])
            hNa[j] += dhNa(Vprev[j], hNa[j])
            mCaT[j] += dmCaT(Vprev[j], mCaT[j])
            hCaT[j] += dhCaT(Vprev[j], hCaT[j])
            mCaS[j] += dmCaS(Vprev[j], mCaS[j])
            hCaS[j] += dhCaS(Vprev[j], hCaS[j])
            mA[j] += dmA(Vprev[j], mA[j])
            hA[j] += dhA(Vprev[j], hA[j])
            mKd[j] += dmKd(Vprev[j], mKd[j])
            mH[j] += dmH(Vprev[j], mH[j])

            # Computing new reference values of gCaS and gA if DICs have changed
            if !((gs_th[z, j] == gs_th[z-1, j]) & (gu_th[z, j] == gu_th[z-1, j]))
                (gCaS_r[j], gA_r[j]) = DICs_gmax_neuromodCaSA(gNa[j], gCaT[j], gKd[j], gKCa[j], gH[j],
                                                    gleak[j], gs_th[z, j], gu_th[z, j], Vth[j], tmKCa[j])
            end

            # Computing control signals
            uCaS, vCaS, eCaS = uve(gCaS[j], gCaS_r[j], zCaS[j], Kp, Ki, u_maxCaS)
            uA, vA, eA = uve(gA[j], gA_r[j], zA[j], Kp, Ki, u_maxA)

            # Controller ODEs
            gCaSi[j] += dgi(gCaS[j], gCaSi[j], α, β, uCaS)
            gCaS[j] += dg(gCaS[j], gCaSi[j], α)
            zCaS[j] += dz(eCaS, Kt, uCaS, vCaS)
            gAi[j] += dgi(gA[j], gAi[j], α, β, uA)
            gA[j] += dg(gA[j], gAi[j], α)
            zA[j] += dz(eA, Kt, uA, vA)
        end

        # Coupling ODEs
        Vcur[1] += (dt) * (1/C) * (-gsyn21*s21*(Vprev[1]-Vsyn) - gsyn51*s51*(Vprev[1]-Vsyn) - gsyn81*s81*(Vprev[1]-Vsyn))
        s21 += ds(Vprev[2], s21)
        s51 += ds(Vprev[5], s51)
        s81 += ds(Vprev[8], s81)

        Vcur[2] += (dt) * (1/C) * (-gsyn12*s12*(Vprev[2]-Vsyn) - gsyn62*s62*(Vprev[2]-Vsyn) - gsyn72*s72*(Vprev[2]-Vsyn))
        s12 += ds(Vprev[1], s12)
        s62 += ds(Vprev[6], s62)
        s72 += ds(Vprev[7], s72)

        Vcur[3] += (dt) * (1/C) * (-gsyn43*s43*(Vprev[3]-Vsyn) - gsyn53*s53*(Vprev[3]-Vsyn) - gsyn73*s73*(Vprev[3]-Vsyn))
        s43 += ds(Vprev[4], s43)
        s53 += ds(Vprev[5], s53)
        s73 += ds(Vprev[7], s73)

        Vcur[4] += (dt) * (1/C) * (-gsyn34*s34*(Vprev[4]-Vsyn) - gsyn64*s64*(Vprev[4]-Vsyn) - gsyn84*s84*(Vprev[4]-Vsyn))
        s34 += ds(Vprev[3], s34)
        s64 += ds(Vprev[6], s64)
        s84 += ds(Vprev[8], s84)

        Vcur[5] += (dt) * (1/C) * (-gsyn65*s65*(Vprev[5]-Vsyn))
        s65 += ds(Vprev[6], s65)

        Vcur[6] += (dt) * (1/C) * (-gsyn56*s56*(Vprev[6]-Vsyn))
        s56 += ds(Vprev[5], s56)

        Vcur[7] += (dt) * (1/C) * (-gsyn87*s87*(Vprev[7]-Vsyn))
        s87 += ds(Vprev[8], s87)

        Vcur[8] += (dt) * (1/C) * (-gsyn78*s78*(Vprev[8]-Vsyn))
        s78 += ds(Vprev[7], s78)

        Vprev = copy(Vcur)
        V_sol[z, :] = copy(Vcur)
    end

    return V_sol
end

## Stimulation function
heaviside(t) = (1 + sign(t)) / 2
pulse(t, ti, tf) = heaviside(t-ti) - heaviside(t-tf)
