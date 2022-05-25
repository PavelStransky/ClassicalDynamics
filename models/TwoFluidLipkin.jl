using Random
using IntervalRootFinding

""" Classical energy of the Two-Fluid Lipkin model """
function Energy(x, parameters)
    p1, p2, q1, q2 = x
    ξ, η1, η2 = parameters

    Σ1 = p1*p1 + q1*q1 
    Σ2 = p2*p2 + q2*q2
    
    s1 = 2 - Σ1
    s2 = 2 - Σ2

    if s1 < 0.0 || s2 < 0
        return NaN
#        print("s negative!")
#        s2 = 0.0
    end

    s1 = sqrt(s1)
    s2 = sqrt(s2)

    return 0.5 * ξ * (Σ1 + Σ2) - 0.125 * (1 - ξ) * (2 * q1 * s1 + 2 * q2 * s2 + η1 * Σ1 + η2 * Σ2)^2
end

""" Filling initial conditions (at a given energy calculates randomly all missing points) """
function InitialConditions(x0, e, parameters, coordinate)
    x = deepcopy(x0)
    x[coordinate] = 0

    p1, p2, q1, q2 = x[1:4]
    ξ, η1, η2 = parameters

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
#        print("s negative!")
#        s2 = 0.0
    end

    s1 = sqrt(s1)
    s2 = sqrt(s2)

    Φ = reshape(x[5:end], 4, 4)     # Tangent dynamics

    ξ, η1, η2 = parameters.modelParameters

    Ξ = 0.5 * (1 - ξ)

    A1 = -q1 * q1 / s1 + s1 + η1 * q1
    A2 = -q2 * q2 / s2 + s2 + η2 * q2
    B = 2 * q1 * s1 + 2 * q2 * s2 + η1 * Σ1 + η2 * Σ2

    c1 = η1 - q1 / s1
    C1 = p1 * c1

    c2 = η2 - q2 / s2
    C2 = p2 * c2

    r1 = q1 / s1
    r2 = q2 / s2

    dx[1] = -ξ * q1 + Ξ * A1 * B
    dx[2] = -ξ * q2 + Ξ * A2 * B
    dx[3] = ξ * p1 - Ξ * C1 * B
    dx[4] = ξ * p2 - Ξ * C2 * B

    Hq1p1 = -Ξ * (2 * A1 * C1 - p1 / s1 * (r1^2 + 1) * B)
    Hq1p2 = -2 * Ξ * A1 * C2
    Hq1q1 = ξ - Ξ * (2 * A1 * A1 + (-r1^3 - 3 * r1 + η1) * B)
    Hq1q2 = -2 * Ξ * A1 * A2

    Hq2p1 = -2 * Ξ * A2 * C1
    Hq2p2 = -Ξ * (2 * A2 * C2 - p2 / s2 * (r2^2 + 1) * B)
    Hq2q2 = ξ - Ξ * (2 * A2 * A2 + (-r2^3 - 3 * r2 + η2) * B)

    Hp1p1 = ξ - Ξ * (2 * C1 * C1 + ((-p1 * p1 * q1) / s1^3 + c1) * B)
    Hp1p2 = -2 * Ξ * C1 * C2

    Hp2p2 = ξ - Ξ * (2 * C2 * C2 + ((-p2 * p2 * q2) / s2^3 + c2) * B)
    
    j = [-Hq1p1 -Hq1p2 -Hq1q1 -Hq1q2;
         -Hq2p1 -Hq2p2 -Hq1q2 -Hq2q2;
         Hp1p1 Hp1p2 Hq1p1 Hq2p1;
         Hp1p2 Hp2p2 Hq1p2 Hq2p2]

    dx[5:end] = j * Φ
end
