using Random
using Logging
using IntervalRootFinding

# x = (p, q)

""" Energy of the Bose-Hubbard model """
function Energy(x, parameters)
    L, J, U = parameters

    energy = 0.0
    for i = 1:L
        j = i + L
        k = i == L ? 1 : i + 1
        l = i == L ? L + 1 : j + 1

        energy += -J * (x[i] * x[k] + x[j] * x[l])
        energy += U / 4 * (x[i] * x[i] + x[j] * x[j]) ^ 2
    end

    return energy
end

""" True if given trajectory is within the kinematically accessible domain """
function CheckDomain(x, parameters, t)
    return false
end

function InitialCondition(energy, parameters, error; maxInitialConditions=1000000)
    L, J, U = parameters

    for i = 1:maxInitialConditions
        x = randn(2 * L)
        x .*= sqrt.(2 / sum(x .* x))

        e = Energy(x, parameters)

        if abs(energy - e) < error
            @info "Initial condition after $(i) attempts found: $(x)"
            return x
        end
    end

    return nothing
end


function EquationOfMotion!(dx, x, parameters, t)
    L, J, U = parameters.modelParameters

    # Trajectories
    for i = 1:L
        a = U * (x[i] * x[i] + x[i + L] * x[i + L])

        qm1 = i == 1 ? x[2 * L] : x[i + L - 1]
        qp1 = i == L ? x[L + 1] : x[i + L + 1]        
        dx[i] = J * (qm1 + qp1) - a * x[i + L]

        pm1 = i == 1 ? x[L] : x[i - 1]
        pp1 = i == L ? x[1] : x[i + 1]
        dx[i + L] = -J * (pm1 + pp1) + a * x[i]
    end

    Φ = reshape(x[(2 * L + 1):end], 2 * L, 2 * L)     # Tangent dynamics
    j = zeros(Float64, 2 * L, 2 * L)

    for i = 1:(2 * L)
        for k = 1:(2 * L)
            i1 = i
            i2 = (k + L - 1) % (2 * L) + 1
            sgn = i2 > L ? -1 : 1

            v = 0
            if i1 == i2
                v = U * (3 * x[i] * x[i] + x[k] * x[k])
            elseif (i1 <= L && i2 <= L && (abs(i1 - i2) == 1 || abs(i1 - i2) == L - 1)) || (i1 > L && i2 > L && (abs(i1 - i2) == 1 || abs(i1 - i2) == L - 1))
                v = -J
            elseif i == k
                v = 2 * U * x[i1] * x[i2]
            end

            j[k, i] = sgn * v
        end
    end

    # println(x[1:(2 * L)])
    # println("j = ", j)
    # println(Energy(x, parameters.modelParameters), " at t = ", t)

    # j = zeros(Float64, 2 * L, 2 * L)
    dx[(2 * L + 1):end] = j * Φ
end
