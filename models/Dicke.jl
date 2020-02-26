using Random

""" Energy of the Dicke model """
function Energy(x, parameters)
    P, p, Q, q = x
    λ, δ, ω, ω₀  = parameters

    return 0.5 * ω * (p*p + q*q) + 0.5 * ω₀ * (P*P + Q*Q) + λ * sqrt(1.0 - 0.25 * (P*P + Q*Q)) * ((1.0 + δ) * P * q - (1.0 - δ) * p * Q)
end

""" Filling initial conditions (at a given energy calculates randomly all missing points) """
function InitialCondition!(x0, e, parameters)
    numMissingValues = length(x0) - length(collect(skipmissing(x0)))

    P, p, Q, q = x0
    λ, δ, ω, ω₀  = parameters

    if numMissingValues == 1 && ismissing(q) && p == 0.0
        s2 = 1.0 - 0.25 * (P*P + Q*Q)
        if s2 < 0
            return false
        end
    
        s = sqrt(s2)
        b = (1.0 + δ) * λ * P * s
    
        discriminant = b*b + 2.0 * ω * (e - 0.5 * ω₀ * (P*P + Q*Q))
        if discriminant < 0
            return false
        end

        sign = rand() >= 0.5 ? 1.0 : -1.0

        q = (sign * sqrt(discriminant) - b) / ω
        x0[4] = q

        return true
    end

    if numMissingValues == 1 && ismissing(Q) && P == 0.0
        a = λ * (1.0 - δ) * p
        b = 0.5 * ω * (p * p + q * q)

        discriminant1 = a * a - (e - b) * (e - 2.0 * ω₀ - b)

        if discriminant1 < 0
            return false
        end

        c = e * ω₀ + a * a - ω₀ * b
        d = ω₀ * ω₀ + a * a

        discriminant2p = 2.0 * (c + a * sqrt(discriminant1)) / d
        discriminant2m = 2.0 * (c - a * sqrt(discriminant1)) / d

        if discriminant2p < 0 && discriminant2m < 0
            return false
        end

        # Not all of the found initial conditions have the desired energy; let's select the correct ones
        ics = []

        x0[3] = sqrt(discriminant2m)
        if isapprox(e, Energy(x0, parameters))
            append!(ics, x0[3])
        end

        x0[3] = -sqrt(discriminant2m)
        if isapprox(e, Energy(x0, parameters))
            append!(ics, x0[3])
        end
    
        x0[3] = sqrt(discriminant2p)
        if isapprox(e, Energy(x0, parameters))
            append!(ics, x0[3])
        end
    
        x0[3] = -sqrt(discriminant2p)
        if isapprox(e, Energy(x0, parameters))
            append!(ics, x0[3])
        end

        println(ics)

        if length(ics) == 0
            x0[3] = missing
            return false
        end        

        x0[3] = rand(ics)
        return true
    end

    return false
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
    λ, δ, ω, ω₀ = parameters

    s = sqrt(s2)
    b = (1.0 + δ) * P * q - (1.0 - δ) * p * Q

    ls_4 = 0.25 * λ * s
    l_4s = 0.25 * λ/ s

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