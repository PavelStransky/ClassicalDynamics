j = 30
maxn = 70
ω = 1.0
ω₀ = 1.0
δ = 1.0

λ = 1.0
λᵪ= sqrt(ω * ω₀) / (1.0 + δ)

G = 2.0 * λ / (sqrt(2.0 * j) * ω)
logG = log(G)

function Index(n, m)
    return n * (2.0 * j + 1) + m + j + 1
end

function LogFactorial(n)
    if n <= 1
        return 0
    else
        return log(n) + LogFactorial(n - 1)
    end
end

function Minus(n)
    if n % 2 == 0
        return 1
    else
        return -1
    end
end

function fu(n1, n2)
    coef = 0.5 * (LogFactorial(n1) + LogFactorial(n2) - G*G)

    sum = 0
    for r = 0:min(n1, n2)
        sum += Minus(r) * exp(logG * (n1 + n2 - 2.0 * r) - LogFactorial(n1 - r) - LogFactorial(n2 - r) - LogFactorial(r))
    end

    return sum
end

function Fuf(np, mp, n, m)
    if m + 1 == mp
        return -0.5 * ω₀ * sqrt(j * (j + 1) - m * (m + 1)) * Minus(n) * fu(np + 1, n + 1)
    elseif m - 1 == mp
        return -0.5 * ω₀ * sqrt(j * (j + 1) - m * (m - 1)) * Minus(np) * fu(np + 1, n + 1)
    end

    return 0.0
end
