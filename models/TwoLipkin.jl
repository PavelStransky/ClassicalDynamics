using Random
using IntervalRootFinding

""" Classical energy of the Two Lipkin (Šindelka) model """
function Energy(x, parameters)
    p1, p2, q1, q2 = x
    ω, λ, δ = parameters

    Σ1 = p1*p1 + q1*q1 
    Σ2 = p2*p2 + q2*q2
    
    s1 = 2 - Σ1
    s2 = 2 - Σ2

    if s1 < 0.0 || s2 < 0
        return NaN
    end

    s1 = sqrt(s1)
    s2 = sqrt(s2)

    return 0.5 * ω * (Σ1 + Σ2) + 0.25 * λ * s1 * s2 * (q1 * q2 + p1 * p2 + δ * (q1 * q2 - p1 * p2))
end

""" Filling initial conditions (at a given energy calculates randomly all missing points) """
function InitialConditions(x0, e, parameters, coordinate)
    x = deepcopy(x0)
    x[coordinate] = 0

    p1, p2, q1, q2 = x[1:4]
    ω, λ, δ = parameters

    result = []
    
    Σ1 = p1*p1 + q1*q1 
    Σ2 = p2*p2 + q2*q2
    
    s1 = 2 - Σ1
    s2 = 2 - Σ2

    if s1 >= 0 && s2 >= 0
        try
            xs = LinRange(-sqrt(2), sqrt(2), 1000)

            if coordinate == 1
                ys = [Energy((x, p2, q1, q2), parameters) for x in xs]
                xs = [x for (x, y) in zip(xs, ys) if !isnan(y)]
                result = roots(x -> Energy((x, p2, q1, q2), parameters) - e, minimum(xs)..maximum(xs))
            
            elseif coordinate == 2
                ys = [Energy((p1, x, q1, q2), parameters) for x in xs]
                xs = [x for (x, y) in zip(xs, ys) if !isnan(y)]
                result = roots(x -> Energy((p1, x, q1, q2), parameters) - e, minimum(xs)..maximum(xs))
            
            elseif coordinate == 3
                ys = [Energy((p1, p2, x, q2), parameters) for x in xs]
                xs = [x for (x, y) in zip(xs, ys) if !isnan(y)]
                result = roots(x -> Energy((p1, p2, x, q2), parameters) - e, minimum(xs)..maximum(xs))

            elseif coordinate == 4
                ys = [Energy((p1, p2, q1, x), parameters) for x in xs]
                xs = [x for (x, y) in zip(xs, ys) if !isnan(y)]
                result = roots(x -> Energy((p1, p2, q1, x), parameters) - e, minimum(xs)..maximum(xs))
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
    p1, p2, q1, q2 = x[1:4]

    Σ1 = p1*p1 + q1*q1 
    Σ2 = p2*p2 + q2*q2
    
    s1 = 2 - Σ1
    s2 = 2 - Σ2

    return s1 < 0 || s2 < 0
end

function EquationOfMotion!(dx, x, parameters, t)
    p1, p2, q1, q2 = x[1:4]             # Equations of motion for trajectories

    Σ1 = p1*p1 + q1*q1 
    Σ2 = p2*p2 + q2*q2
    
    s1 = 2 - Σ1
    s2 = 2 - Σ2

    if s1 < 0.0 || s2 < 0
        return NaN
    end

    s1 = sqrt(s1)
    s2 = sqrt(s2)

    Φ = reshape(x[5:end], 4, 4)     # Tangent dynamics

    ω, λ, δ = parameters.modelParameters

    λ *= 0.25

    c = q1 * q2 + p1 * p2 + δ * (q1 * q2 - p1 * p2)

    δp = 1 + δ
    δm = 1 - δ

    dx[1] = -ω * q1 - λ * s2 * (s1 * δp * q2 - c / s1 * q1)
    dx[2] = -ω * q2 - λ * s1 * (s2 * δp * q1 - c / s2 * q2)
    dx[3] = ω * p1 + λ * s2 * (s1 * δm * p2 - c / s1 * p1)
    dx[4] = ω * p2 + λ * s1 * (s2 * δm * p1 - c / s2 * p2)

    Hq1p1 = -λ * s2 / s1 * (δm * q1 * p2 + δp * q2 * p1 + q1 * p1 * c / s1^2)
    Hq1p2 = -λ * (s2 / s1 * δm * q1 * p1 + s1 / s2 * δp * s2 * p2 - c / (s1 * s2) * q1 * p2)
    Hq1q1 = ω - λ * s2 / s1 * (2 * δp * q1 * q2 + c / s1^2 * q1^2 + c) 
    Hq1q2 = -λ * (-s1 * s2 * δp + s2 / s1 * δp * q1^2 + s1 / s2 * δp * q2^2 - c / (s1 * s2) * q1 * q2)

    Hq2p1 = -λ * (s2 / s1 * δp * q1 * p1 + s1 / s2 * δm * s2 * p2 - c / (s1 * s2) * q2 * p1)
    Hq2p2 = -λ * s1 / s2 * (δp * q1 * p2 + δm * q2 * p1 + q2 * p2 * c / s2^2)
    Hq2q2 = ω - λ * s1 / s2 * (2 * δp * q1 * q2 + c / s2^2 * q2^2 + c) 

    Hp1p1 = ω - λ * s2 / s1 * (2 * δm * p1 * p2 + c / s1^2 * p1^2 + c) 
    Hp1p2 = -λ * (-s1 * s2 * δm + s2 / s1 * δm * p1^2 + s1 / s2 * δm * p2^2 - c / (s1 * s2) * p1 * p2)

    Hp2p2 = ω - λ * s1 / s2 * (2 * δm * p1 * p2 + c / s2^2 * p2^2 + c) 
    
    j = [-Hq1p1 -Hq1p2 -Hq1q1 -Hq1q2;
         -Hq2p1 -Hq2p2 -Hq1q2 -Hq2q2;
         Hp1p1 Hp1p2 Hq1p1 Hq2p1;
         Hp1p2 Hp2p2 Hq1p2 Hq2p2]

    dx[5:end] = j * Φ
end
