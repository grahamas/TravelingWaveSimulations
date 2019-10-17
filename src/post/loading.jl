function binarize_directory(dir, chunks)
    csv = joinpath(dir, "*.csv")
    bin = joinpath(dir, "bin")
    loadtable(glob(csv), output=bin, chunks=chunks)
end

function load_directory(dir, chunks=-1)
    bin_dir = joinpath(dir, "bin")
    if chunks < 1
        if !isdir(bin_dir)
             error("Directory not binarized; specify chunks.")
        end
    else
        binarize_directory(dir, chunks)
    end
    return load(bin_dir)
end
