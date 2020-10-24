#!/home/grahams/.local/bin/julia
#SBATCH --mail-user=remote.compute.results@gmail.com 
#SBATCH --mail-type=END
#SBATCH --ntasks=24 
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=fast 
#SBATCH --time=1-00:00:00

n_procs = 24
data_root = "/scratch/grahams"

using DrWatson
quickactivate("/home/grahams/git/TravelingWaveSimulations")
include(projectdir("_research", "repl", "parallel", "init_n_procs.jl"))
include(projectdir("_research", "repl", "depblock_strength_threshold_sweep.jl"))
