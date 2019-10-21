function binarize_directory(dir, chunks)
    csv = joinpath(dir, "*.csv")
    bin = joinpath(dir, "bin")
    loadtable(glob(csv), output=bin, chunks=chunks)
end

function load_directory(dir, args...)
    if ispath(joinpath(dir, "juliadb_index"))
        return load(dir, args...)
    end
    bin_dir = joinpath(dir, "bin")
    if ispath(joinpath(bin_dir, "juliadb_index")) 
        return load(bin_dir)
    else
        binarize_directory(dir, args...)
        return load(bin_dir)
    end
end
