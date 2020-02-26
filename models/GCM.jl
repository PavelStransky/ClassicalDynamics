using Random

""" Energy of the GCM model """
function Energy(X, parameters)
    x, y, px, py = X
    A, B, C, K  = parameters

    b2 = x * x + y * y
    p2 = px * px + py * py

    return 0.5 * p2 / K + A * b2 + B * x * (x * x - 3.0 * y * y) + C * b2 * b2
end

""" Filling initial conditions (first for section (x, px) through plane y = 0)"""
function InitialCondition!(x0, e, parameters)
    numMissingValues = length(x0) - length(collect(skipmissing(x0)))

    x, y, px, py = x0
    A, B, C, K  = parameters

    if numMissingValues == 1 && ismissing(py)
        d = 2.0 * K * (e - Energy([x, y, px, 0], parameters))
        if d < 0
            return false
        end

        sign = rand() >= 0.5 ? 1.0 : -1.0
        py = sign * sqrt(d)
        x0[4] = py

        return true
    end

    return false
end

""" True if given trajectory is within the kinematically accessible domain """
function CheckDomain(x, parameters, t)
    return false
end

function EquationOfMotion!(dx, X, parameters, t)
    x, y, px, py = X[1:4]             # Equations of motion for trajectories
    Φ = reshape(X[5:end], 4, 4)     # Tangent dynamics

    A, B, C, K = parameters[1:4]

    b2 = x * x + y * y    
    Hx = 2.0 * A * x + 3.0 * B * (x * x - y * y) + 4.0 * C * x * b2
    Hy = 2.0 * y * (A - 3.0 * B * x + 2.0 * C * b2)

    Hxx = 2.0 * (A + 3.0 * B * x + 2.0 * C * (3.0 * x * x + y * y))
    Hxy = y * (-6.0 * B + 8.0 * C * x)
    Hyy = 2.0 * (A - 3.0 * B * x + 2.0 * C * (x * x + 3.0 * y * y))

    dx[1] = px / K
    dx[2] = py / K
    dx[3] = -Hx
    dx[4] = -Hy

    j = [0 0 1.0/K 0;
         0 0 0 1.0/K;
         -Hxx -Hxy 0 0;
         -Hxy -Hyy 0 0]

    dx[5:end] = j * Φ
end