module RichardsonExtrapolationUQ

using Statistics
using FinEtools
using FinEtools.AlgoBaseModule: richextrapol

function richextrapol_uq(solutions, elementsizes; W = 2/3)
    @assert length(solutions) > 3
    @assert length(elementsizes) > 3
    @assert length(elementsizes)  == length(solutions)
    # Four possible combinations of results
    c = [[1, 2, 3], [1, 3, 4], [1, 2, 4], [2, 3, 4]]

    results = []
    for j in 1:length(solutions)-3
        ess = elementsizes[j:j+3]
        qs = solutions[j:j+3]

        extrapolations = []
        for j in 1:length(c)
            e = (Inf * sign(qs[end]), 0.0, 0.0, Inf)
            try
                e = richextrapol(qs[c[j]], ess[c[j]])
                if e[4] > minimum(abs.(solutions)) / 1.0e-6
                    error("Richardson extrapolation failed: large residual ($e[4])")
                end 
            catch
            end
            println("extrapolation $(qs[c[j]]) $(e[1])")
            println("convergence rate $(e[2])")
            push!(extrapolations, (solnestim = e[1], beta = e[2], c = e[3], maxresidual = e[4]))
        end
        extrsols = [e.solnestim for e in extrapolations]
        q_m = median(extrsols)
        q_m_ad = 1.4826 * median(abs.(extrsols .- q_m))
        q_star = W * extrapolations[end].solnestim + (1 - W) * q_m
        beta_m = median([e.beta for e in extrapolations])
        beta_m_ad = median(abs.([e.beta for e in extrapolations] .- beta_m))
        beta_star = W * extrapolations[end].beta + (1 - W) * beta_m
        push!(results, (estim = q_star, estim_ad_x_2 = 2*q_m_ad, beta = beta_star, beta_ad_x_2 = 2*beta_m_ad, elementsize = ess[end-1], extrapolations = extrapolations)) 
    end 
    return results
end 

function gci_uq(solutions, elementsizes)
    @assert length(solutions) > 3
    @assert length(elementsizes) > 3
    @assert length(elementsizes)  == length(solutions)
    
    results = []
    for j in 1:length(solutions)-2
        ess = elementsizes[j:j+2]
        qs = solutions[j:j+2]

        e = richextrapol(qs, ess)
        if e[4] > minimum(abs.(solutions)) / 1.0e-6
            error("Richardson extrapolation failed: large residual ($e[4])")
        end 
        w_1 = qs[end]
        w_2 = qs[end-1]
        gamma_1 = ess[end]
        gamma_2 = ess[end-1]
        p = e[2]
        println("convergence rate $(p)")
        
        e_a12 = (w_1 - w_2) / w_1 / ((gamma_2 / gamma_1)^p - 1)
        push!(results, (w_1 = w_1, gci = 1.25 * e_a12, beta = p, elementsize = gamma_1)) 
    end 
    return results
end 

# Effective Convergence Checks for Verifying Finite Element Stresses at
# Three-Dimensional Stress Concentrations J. R. Beisheim Development
# Department, ANSYS, Inc., Canonsburg, PA 15317 e-mail: jeff.beisheim@ansys.com
# G. B. Sinclair Department of Mechanical Engineering, Louisiana State
# University, Baton Rouge, LA 70803 e-mail: sinclair@lsu.edu P. J. Roache
# Consultant, 1215 Apache Drive, Socorro, NM 87801 e-mail: hermosa@sdc.org 
function bsr_uq(solutions, elementsizes)
    @assert length(solutions) > 3
    @assert length(elementsizes) > 3
    @assert length(elementsizes)  == length(solutions)
    
    results = []
    for j in 1:length(solutions)-2
        ess = elementsizes[j:j+2]
        qs = solutions[j:j+2]

        e = richextrapol(qs, ess)
        if e[4] > minimum(abs.(solutions)) / 1.0e-6
            error("Richardson extrapolation failed: large residual ($e[4])")
        end 
        @show lambda_m1 = ess[1] / ess[2]
        @show lambda_m = ess[2] / ess[3]
        Delta_q_m1 = qs[2] - qs[1]
        Delta_q_m = qs[3] - qs[2]
        @show chat_m = log(Delta_q_m1 / Delta_q_m) / log(lambda_m)
        println("convergence rate $(chat_m) vs. $(e[2])")
        epshat_m = abs(Delta_q_m) / (Delta_q_m1 / Delta_q_m - 1) / abs(qs[3])
        println("percentage error $(epshat_m)")
        
        push!(results, (q_m = qs[3], epshat_m = epshat_m, chat_m = chat_m, elementsize = ess[3])) 
    end 
    return results
end 

end # module