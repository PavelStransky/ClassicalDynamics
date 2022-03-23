using Random
using IntervalRootFinding

# px = p 2
# py = P 1
# x = q  4
# y = Q  3

""" Energy of the Dicke model """
function Energy(x, parameters)
    P, p, Q, q = x
    A, B, C = parameters

    Σ2 = p*p + P*P + q*q + Q*Q
    s2 = 2.0 - Σ2

    if s2 < 0.0
        return NaN
#        print("s negative!")
#        s2 = 0.0
    end

    return A * Σ2 + B * ((P*P + p*p) * s2 + (P*q - p*Q) ^ 2) + C * P * sqrt(s2)
end

""" Filling initial conditions (at a given energy calculates randomly all missing points) """
function InitialConditions(x0, e, parameters, coordinate)
    x = deepcopy(x0)
    x[coordinate] = 0

    P, p, Q, q = x
    A, B, C  = parameters

    result = []
    
    s2 = 2.0 - (P*P + p*p + Q*Q + q*q)
    if s2 >= 0
        try
            xs = LinRange(-sqrt(2), sqrt(2), 1000)

            if coordinate == 1
                ys = [Energy((x, p, Q, q), parameters) for x in xs]
                xs = [x for (x, y) in zip(xs, ys) if !isnan(y)]
                result = roots(x -> Energy((x, p, Q, q), parameters) - e, minimum(xs)..maximum(xs))
            
            elseif coordinate == 2
                ys = [Energy((P, x, Q, q), parameters) for x in xs]
                xs = [x for (x, y) in zip(xs, ys) if !isnan(y)]
                result = roots(x -> Energy((P, x, Q, q), parameters) - e, minimum(xs)..maximum(xs))
            
            elseif coordinate == 3
                ys = [Energy((P, p, x, q), parameters) for x in xs]
                xs = [x for (x, y) in zip(xs, ys) if !isnan(y)]
                result = roots(x -> Energy((P, p, x, q), parameters) - e, minimum(xs)..maximum(xs))

            elseif coordinate == 4
                ys = [Energy((P, p, Q, x), parameters) for x in xs]
                xs = [x for (x, y) in zip(xs, ys) if !isnan(y)]
                result = roots(x -> Energy((P, p, Q, x), parameters) - e, minimum(xs)..maximum(xs))
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
    s2 = 2.0 - (p*p + P*P + q*q + Q*Q)
    return s2 < 0
end

function EquationOfMotion!(dx, x, parameters, t)
    P, p, Q, q = x[1:4]             # Equations of motion for trajectories

    s2 = 2.0 - (p*p + P*P + q*q + Q*Q)
    if s2 < 0                       
        s2 = 0.0                    # With isoutofdomain=CheckDomain the calculation shouldn't enter here, but in enters anyway
        print("s negative!") 
    end
    
    Φ = reshape(x[5:end], 4, 4)     # Tangent dynamics
    A, B, C = parameters.modelParameters

    s = sqrt(s2)
    a = p*p + P*P
    b = P*Q + p*q

    A2 = 2.0 * A
    B2 = 2.0 * B
    Cs = C / s
    Cs3 = Cs / s2

    dx[1] = -A2 * Q + B2 * P * b + Cs * P * Q
    dx[2] = -A2 * q + B2 * p * b + Cs * P * q
    dx[3] = A2 * P - B2 * (2.0 * P * (a - 1.0) + Q * b) + Cs * (s2 - P*P)
    dx[4] = A2 * p - B2 * (2.0 * p * (a - 1.0) + q * b) - Cs * p * P

    HPQ = -B2 * (p * q + 2.0 * P * Q) - Cs3 * Q * (s2 + P*P)
    HpQ = -B2 * P * q - Cs3 * p * P * Q
    HQQ = A2 - B2 * P*P - Cs3 * P * (s2 + Q*Q)
    HQq = -B2 * p * P - Cs3 * P * q * Q

    HPq = -B2 * p * Q - Cs3 * q * (s2 + P*P)
    Hpq = -B2 * (2.0 * p * q + P * Q) - Cs3 * p * P * q
    Hqq = A2 - B2 * p*p - Cs3 * P * (s2 + q*q)

    HPP = A2 - B2 * (6.0 * P*P + 2.0 * p*p + Q*Q - 2.0) - Cs3 * P * (P*P - 3.0 * a + 3.0 * (s2 + a))
    HPp = -B2 * (4.0 * p * P + q * Q) - Cs3 * p * (s2 + P*P)

    Hpp = A2 - B2 * (2.0 * P*P + 6.0 * p*p + q*q - 2.0) - Cs3 * P * (s2 + p*p)
    
    j = [-HPQ -HpQ -HQQ -HQq;
         -HPq -Hpq -HQq -Hqq;
         HPP HPp HPQ HPq;
         HPp Hpp HpQ Hpq]

    dx[5:end] = j * Φ
end
