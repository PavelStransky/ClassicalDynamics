# run_laptop_number_conserving.ps1
#
# Runs the LAPTOP part of number-conserving-BH-paper-TODO.md, stage by stage.  The split between
# this laptop (i9-13900HX, 24 cores, 32 GB), IPNP36 (24 cores, 640 GB) and the Chimera cluster, the
# cost of every stage and what each result feeds is in number-conserving-BH-compute-plan.md.
#
# By that plan the laptop runs checks, k1, transient, attractors and analysis.  The quantum stages
# (reference, cuts, plane, eta0, variants, dsff, modes, steady, trajectories, sparse) run on IPNP36
# at larger N through run_ipnp36_number_conserving.sh; they stay here, sized for 32 GB, as a
# fallback.  The classical Lyapunov sweeps go to Chimera through BHMapNumberConservingChimeraSubmit.sh.
#
# Usage (PowerShell, from the repository root):
#   .\run_laptop_number_conserving.ps1 -Stage checks
#   .\run_laptop_number_conserving.ps1 -Stage reference, cuts, transient
#   .\run_laptop_number_conserving.ps1 -Stage list
#
# Every stage is resumable: the Julia drivers skip what is already on disk, so an interrupted stage
# is simply started again.  Output goes to ~/results/bh/number-conserving/..., a log of every stage
# to ~/results/bh/number-conserving/logs/<stage>.log.
#
# Stages (wall-clock estimates on this laptop; run the long ones overnight):
#   checks        consistency checks of every new piece of code, sparse vs dense         ~ 10 min
#   k1            K1   master equation vs classical flow, N = 10, 20, 30, 40              ~ 30 min
#   reference     C4 C6 K3 K8  spectra at the four reference points, N = 6 ... 16          ~  2 h
#   cuts          C1 C6 K8  eta = 3 and g = -20 cuts, N = 8 ... 14                        ~  2 h
#   plane         C1   (g, eta) plane at N = 12, clean sector, 14641 points              ~  5 h
#   eta0          C2   eta = 0 line with reflection parity, N = 8 ... 16                  ~  3 h
#   variants      C9 C11 K12  kappa = 0, Gamma_- = Gamma_+/2, kappa = 0.1 and 0.5          ~  2 h
#   dsff          C10  41 spectra around the chaotic point, N = 12                        ~  5 min
#   transient     C13  transient chaos, eta = 0 and eta = 3                               ~  1 h
#   attractors    C14 C15 C16 C17 K9 K10 K11                                              ~  2 h
#   modes         C5   observable weights with eigenvectors, N = 8 ... 14                ~  2 h
#   steady        C7 A10  steady-state Husimi function vs classical measure, N = 8 ... 24  ~ 30 min
#   trajectories  C8   quantum trajectories at the multistable point, N = 10 ... 70        ~  2 h
#   sparse        C3 C6  slow region and gap at N = 18, 20, 24 (sparse slicing)           ~ 10 h
#   analysis      C4 K4 K7 C10 A1 and every figure (after the stages above)              ~ 30 min

param(
    [string[]]$Stage = @("list"),
    [int]$Workers = 20,          # Distributed workers for the dense Liouvillian map
    [int]$Threads = 20,          # Julia threads for the multithreaded scripts
    [int]$SliceProcesses = 6     # parallel processes for the sparse slices (about 3 GB each at N = 24)
)

$ErrorActionPreference = "Continue"
$Repository = $PSScriptRoot
$Logs = Join-Path $HOME "results\bh\number-conserving\logs"
New-Item -ItemType Directory -Force $Logs | Out-Null
Set-Location $Repository

function Run([string]$Log, [string[]]$Arguments) {
    $file = Join-Path $Logs "$Log.log"
    "==== $(Get-Date -Format s)  julia $($Arguments -join ' ')" | Out-File -Append -Encoding utf8 $file
    & julia @Arguments 2>&1 | ForEach-Object { "$_" } | Tee-Object -FilePath $file -Append
}

function Python([string]$Log, [string[]]$Arguments) {
    $file = Join-Path $Logs "$Log.log"
    "==== $(Get-Date -Format s)  python $($Arguments -join ' ')" | Out-File -Append -Encoding utf8 $file
    & python.exe @Arguments 2>&1 | ForEach-Object { "$_" } | Tee-Object -FilePath $file -Append
}

# Sparse slicing of one run: plan (unless present), slices in parallel processes, merge, and up to
# three refinement rounds until the coverage is complete.
function SparseRun([string]$Name, [string[]]$PlanArguments) {
    $directory = Join-Path $HOME "results\bh\number-conserving\quantum\3\sparse\$Name"
    if (-not (Test-Path (Join-Path $directory "plan.txt"))) {
        Run "sparse" (@("BHNumberConservingLiouvillianSparse.jl", "plan", $Name) + $PlanArguments)
    }
    for ($round = 0; $round -le 3; $round++) {
        $processes = @()
        for ($part = 0; $part -lt $SliceProcesses; $part++) {
            $out = Join-Path $Logs "sparse_$($Name)_part$part.log"
            $processes += Start-Process -FilePath "julia" -PassThru -NoNewWindow `
                -RedirectStandardOutput $out -RedirectStandardError "$out.err" `
                -ArgumentList @("BHNumberConservingLiouvillianSparse.jl", "slices", $Name,
                                "--part", "$part", "--parts", "$SliceProcesses")
        }
        $processes | Wait-Process
        $report = Run "sparse" @("BHNumberConservingLiouvillianSparse.jl", "merge", $Name)
        if (-not ($report -match "INCOMPLETE|missing slices")) { break }
        Run "sparse" @("BHNumberConservingLiouvillianSparse.jl", "refine", $Name) | Out-Null
    }
    Run "sparse" @("BHNumberConservingLiouvillianSparse.jl", "statistics", $Name)
}

$env:CD_NO_PLOTS = "true"

foreach ($s in $Stage) {
    switch ($s) {
        "list" {
            Get-Content $PSCommandPath | Select-String -Pattern "^#   " | ForEach-Object { $_.Line.Substring(1) }
        }
        "checks" {
            # (no `julia -e "..."` here: Windows PowerShell strips the inner quotes of such arguments)
            Run "checks" @("BHNumberConservingChecks.jl")
            Run "checks" @("BHNumberConservingLiouvillianSparse.jl", "validate", "12", "--w", "16")
        }
        "k1" {
            Run "k1" @("-t", "4", "BHNumberConservingQuantumDynamics.jl", "consistency", "--N", "10,20,30,40")
        }
        "reference" {
            Run "reference" @("BHNumberConservingLiouvillianMap.jl", "reference", "--workers", "$Workers")
        }
        "cuts" {
            Run "cuts" @("BHNumberConservingLiouvillianMap.jl", "cut-eta3", "--workers", "$Workers")
            Run "cuts" @("BHNumberConservingLiouvillianMap.jl", "cut-g20", "--workers", "$Workers")
        }
        "plane" {
            Run "plane" @("BHNumberConservingLiouvillianMap.jl", "plane", "--workers", "$Workers")
        }
        "eta0" {
            Run "eta0" @("BHNumberConservingLiouvillianMap.jl", "eta0", "--workers", "$Workers")
        }
        "variants" {
            foreach ($task in @("kappa0", "asymmetric", "kappa")) {
                Run "variants" @("BHNumberConservingLiouvillianMap.jl", $task, "--workers", "$Workers")
            }
        }
        "dsff" {
            Run "dsff" @("BHNumberConservingLiouvillianMap.jl", "dsff", "--workers", "$Workers")
        }
        "transient" {
            Run "transient" @("-t", "$Threads", "BHNumberConservingTransient.jl", "eta0")
            Run "transient" @("-t", "$Threads", "BHNumberConservingTransient.jl", "eta3")
        }
        "attractors" {
            foreach ($command in @("basins", "uncertainty", "bifurcation", "correlations", "simplex", "robustness")) {
                Run "attractors" @("-t", "$Threads", "BHNumberConservingAttractors.jl", $command)
            }
        }
        "modes" {
            Run "modes" @("BHNumberConservingLiouvillianModes.jl", "--N", "8,10,12,14", "--blas", "8")
        }
        "steady" {
            Run "steady" @("-t", "$Threads", "BHNumberConservingSteadyState.jl", "--N", "8,10,12,14,16,20,24")
        }
        "trajectories" {
            Run "trajectories" @("-t", "$Threads", "BHNumberConservingQuantumDynamics.jl", "trajectories",
                                 "--N", "10,20,30,40,50,60,70", "--trajectories", "20", "--time", "2000")
        }
        "sparse" {
            # C3: slow region of the clean sector at the chaotic point (statistics up to w = 16)
            foreach ($N in @(18, 20, 24)) { SparseRun "c3-N$N" @("$N", "1", "21") }
            # C6: gap at the chaotic and at the multistable point, both sectors
            foreach ($point in @(@("-20", "3", "g-20e3"), @("-20", "1", "g-20e1"))) {
                foreach ($N in @(18, 20, 24)) {
                    foreach ($m in @("0", "1")) {
                        SparseRun "c6-$($point[2])-N$N-m$m" @("$N", $m, "3", "--g", $point[0], "--eta", $point[1], "--howmany", "60")
                    }
                }
            }
        }
        "analysis" {
            foreach ($command in @("tables", "unfolding", "edge", "dsff", "classes")) {
                Run "analysis" @("-t", "8", "BHNumberConservingLiouvillianAnalysis.jl", $command, "--N", "12,14,16")
            }
            $figures = Join-Path $HOME "results\bh\number-conserving\figures"
            New-Item -ItemType Directory -Force $figures | Out-Null
            Push-Location $figures
            foreach ($command in @("plane", "cuts", "gap", "eta0")) {
                Python "analysis" @((Join-Path $Repository "analyse_liouvillian_number_conserving.py"), $command, "--no-show")
            }
            foreach ($command in @("transient", "bifurcation", "uncertainty", "simplex", "correlations")) {
                Python "analysis" @((Join-Path $Repository "analyse_attractors_number_conserving.py"), $command, "--no-show")
            }
            Pop-Location
        }
        default { Write-Warning "unknown stage '$s' - run with -Stage list" }
    }
}
