using Random
using IntervalRootFinding

""" Energy of the Dicke model """
function Energy(x, parameters)
    P, p, Q, q = x
    λ, δ, ω, ω₀  = parameters

    s2 = 1.0 - 0.25 * (P*P - Q*Q)
    if s2 < 0.0
        return NaN
    end

    return 0.5 * ω * (p*p + q*q) + 0.5 * ω₀ * (P*P + Q*Q) + λ * sqrt(s2) * ((1.0 + δ) * P * q - (1.0 - δ) * p * Q)
end

""" Filling initial conditions (at a given energy calculates randomly all missing points) """
function InitialConditions(x0, energy, parameters, coordinate)
    x = deepcopy(x0)
    x[coordinate] = 0

    P, p, Q, q = x[1:4]
    λ, δ, ω, ω₀  = parameters

    result = []

    s2 = 1.0 - 0.25 * (P*P + Q*Q)

    if s2 >= 0
        try
            xs = []

            if coordinate == 1 || coordinate == 3 || energy / ω < 0.25
                xs = LinRange(-2, 2, 1000)
            else
                xs = LinRange(-4 * sqrt(energy / ω), 4 * sqrt(energy / ω), 1000)
            end

            if coordinate == 1
                ys = [Energy((x, p, Q, q), parameters) for x in xs]
                xs = [x for (x, y) in zip(xs, ys) if !isnan(y)]
                result = roots(x -> Energy((x, p, Q, q), parameters) - energy, minimum(xs)..maximum(xs))
            
            elseif coordinate == 2
                ys = [Energy((P, x, Q, q), parameters) for x in xs]
                xs = [x for (x, y) in zip(xs, ys) if !isnan(y)]
                result = roots(x -> Energy((P, x, Q, q), parameters) - energy, minimum(xs)..maximum(xs))
            
            elseif coordinate == 3
                ys = [Energy((P, p, x, q), parameters) for x in xs]
                xs = [x for (x, y) in zip(xs, ys) if !isnan(y)]
                result = roots(x -> Energy((P, p, x, q), parameters) - energy, minimum(xs)..maximum(xs))

            elseif coordinate == 4
                ys = [Energy((P, p, Q, x), parameters) for x in xs]
                xs = [x for (x, y) in zip(xs, ys) if !isnan(y)]
                result = roots(x -> Energy((P, p, Q, x), parameters) - energy, minimum(xs)..maximum(xs))
            end

        catch e
        end
    end

    r = zeros(Float64, 0)
    for i = 1:length(result)
        m = mid(result[i].interval)
        if !isnan(m)
            append!(r, m)
        end
    end

    return sort(r, rev=true)        
end

""" True if given trajectory is within the kinematically accessible domain """
function CheckDomain(x, parameters, t)
    P, p, Q, q = x[1:4]
    s2 = 1.0 - 0.25 * (P*P + Q*Q)
    return s2 < 0
end

function EquationOfMotion!(dx, x, parameters, t)
    P, p, Q, q = x[1:4]             # Equations of motion for trajectories

    s2 = 1.0 - 0.25 * (P*P + Q*Q)
    if s2 < 0                       
        s2 = 0                      # With isoutofdomain=CheckDomain the calculation shouldn't enter here, but in enters anyway
        @error("s negative!") 
    end
    
    Φ = reshape(x[5:end], 4, 4)     # Tangent dynamics
    λ, δ, ω, ω₀ = parameters.modelParameters

    s = sqrt(s2)
    b = (1.0 + δ) * P * q - (1.0 - δ) * p * Q

    ls_4 = 0.25 * λ * s
    l_4s = 0.25 * λ / s

    Q2_4s2 = 0.25 * Q * Q / s2    
    P2_4s2 = 0.25 * P * P / s2

    a = l_4s * P * Q

    HPQ = -l_4s * (0.25 * b * P * Q / s2 + (1.0 + δ) * Q * q - (1.0 - δ) * P * p)
    HpQ = λ * s * (1.0 - δ) * (Q2_4s2 - 1.0)
    HQQ = ω₀ - l_4s * b * (Q2_4s2 + 1.0) + 2.0 * l_4s * (1.0 - δ) * p * Q
    HQq = -a * (1.0 + δ)

    HPq = -λ * s * (1.0 + δ) * (P2_4s2 - 1.0)
    Hpq = 0.0   
    Hqq = ω

    HPP = ω₀ - l_4s * b * (P2_4s2 + 1.0) - 2.0 * l_4s * (1.0 + δ) * P * q
    HPp = a * (1.0 - δ)

    Hpp = ω

    dx[1] = -ω₀ * Q + λ * (0.25 * b * Q / s + (1.0 - δ) * p * s)
    dx[2] = -ω * q - λ * (1.0 + δ) * P * s
    dx[3] = ω₀ * P + λ * (-0.25 * b * P / s + (1.0 + δ) * q * s)
    dx[4] = ω * p - λ * (1.0 - δ) * Q * s

    j = [-HPQ -HpQ -HQQ -HQq;
         -HPq -Hpq -HQq -Hqq;
         HPP HPp HPQ HPq;
         HPp Hpp HpQ Hpq]

    dx[5:end] = j * Φ
end