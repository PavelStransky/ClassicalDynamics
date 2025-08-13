using Random
using IntervalRootFinding

# x = (p, q)

function S2(x)
    return 2.0 - sum(x .* x)
end

""" Energy of the Bose-Hubbard model """
function Energy(x, parameters)
    L, J, U = parameters

    s2 = S2(x[1:(2 * L)])

    if s2 < 0.0
        return NaN
#        print("s negative!")
#        s2 = 0.0
    end

    result = -J * (x[1] + x[L]) + U / 4 * s2 * s2

    for i = 1:L
        result += -J * (x[i] * x[i + 1] + x[i + L] * x[i + L + 1])
        result += U / 4 * (x[i] * x[i] + x[i + L] * x[i + L]) ^ 2
    end

    return result   
end

""" Filling initial conditions (at a given energy calculates randomly all missing points) """
function InitialConditions(x0, e, parameters, coordinate)
    x1 = deepcopy(x0)
    x1[coordinate] = 0

    L, J, U = parameters

    result = []
    
    s2 = S2(x1)
    if s2 >= 0
        try
            xs = LinRange(-sqrt(2), sqrt(2), 1000)

            ys = [x -> begin x1[coordinate] = x; Energy(x1, parameters) end for x in xs]
            xs = [x for (x, y) in zip(xs, ys) if !isnan(y)]
            result = roots(x -> begin x1[coordinate] = x; Energy(x1, parameters) - e end, minimum(xs)..maximum(xs))
        catch e
        end
    end

    r = zeros(Float64, 0)
    for r in result
        m = mid(r.interval)
        if !isnan(m)
            append!(r, m)
        end
    end

    return sort(r, rev=true)
end

""" True if given trajectory is within the kinematically accessible domain """
function CheckDomain(x, parameters, t)
    L, J, U = parameters.modelParameters
    return S2(x[1:(2 * L)]) < 0
end

function EquationOfMotion!(dx, x, parameters, t)
    L, J, U = parameters.modelParameters

    s2 = S2(x[1:(2 * L)])
    if s2 < 0                       
        s2 = 0.0                    # With isoutofdomain=CheckDomain the calculation shouldn't enter here, but in enters anyway
        print("s negative!") 
    end

    s = sqrt(s2)

    # Trajectories
    for i = 1:L
        qm1 = i == 1 ? s : x[i - 1]
        qp1 = i == L ? s : x[i + 1]        
        dx[i] = J * (qm1 + qp1 - x[i] * (x[L + 1] + x[2 * L]) / s) - U * x[i] * ((x[i] * x[i] + x[i + L] * x[i + L])^2 - s2)

        pm1 = i == 1 ? 0 : x[i + L - 1]
        pp1 = i == L ? 0 : x[i + L + 1]
        dx[i + L] = -J * (pm1 + pp1 - x[i + L] * (x[1] + x[L]) / s) + U * x[i + L] * ((x[i] * x[i] + x[i + L] * x[i + L])^2 - s2)
    end

    Φ = reshape(x[(2 * L + 1):end], 2 * L, 2 * L)     # Tangent dynamics
    j = zeros(Float64, 2 * L, 2 * L)

    Q = x[L + 1] + x[2 * L]
    s3 = s * s2

    for i = 1:(2 * L)
        for k = 1:(2 * L)
            i1 = i
            i2 = (k + L) % (2 * L)
            sgn = i2 <= L ? -1 : 1
            m = i1 % L == i2 % L ? 4 : 2

            v = x[i1] * x[i2] * (Q / s3 + U * m)

            if i1 == i2
                v += J * Q / s + U * (x[i1] * x[i1] + x[i2] * x[i2] - s2)

                if i1 == L + 1 || i1 == 2 * L
                    v += 2 * J * x[i1] / s
                end
            elseif (i1 <= L && i2 <= L) || (i1 > L && i2 > L)
                v -= J
            elseif i1 == L + 1 || i1 == 2 * L
                v += J * x[i2] / s
            elseif i2 == L + 1 || i2 == 2 * L
                v += J * x[i1] / s
            end

            j[i, k] = v
        end
    end

    dx[(2 * L + 1):end] = j * Φ
end
