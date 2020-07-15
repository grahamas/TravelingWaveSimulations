function binarize_directory(dir, chunks)
    csv = joinpath(dir, "*.csv")
    bin = joinpath(dir, "bin")
    loadtable(glob(csv), output=bin, chunks=chunks)
end

function has_jdb(dir)
    files = readdir(dir)
    return any([splitext(file)[2] == ".jdb" for file in files])
end

function load_directory(dir, args...)
    if ispath(joinpath(dir, "juliadb_index")) || has_jdb(dir)
        return load(dir, args...)
    end
    bin_dir = joinpath(dir, "bin")
    if ispath(joinpath(bin_dir, "juliadb_index")) || has_jdb(bin_dir)
        return load(bin_dir)
    else
        binarize_directory(dir, args...)
        return load(bin_dir)
    end
end
