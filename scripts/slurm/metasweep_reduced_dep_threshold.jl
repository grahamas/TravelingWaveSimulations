#!/home/grahams/.local/bin/julia
#SBATCH --mail-user=remote.compute.results@gmail.com 
#SBATCH --mail-type=END
#SBATCH --ntasks=8 
#SBATCH --mem-per-cpu=2G
#SBATCH --partition=general 
#SBATCH --time=05:00:00

n_procs = 8
data_root = "/scratch/grahams"

using Pkg
Pkg.activate("/home/grahams/git/TravelingWaveSimulations")
using DrWatson
include(projectdir("_research", "repl", "parallel", "init_n_procs.jl"))
include(projectdir("_research", "repl", "depblock_reduced_threshold_sweep.jl"))
