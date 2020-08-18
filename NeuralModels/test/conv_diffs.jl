
test_ns = 70:10:1000
fft_diffs = zeros(length(test_ns)); naive_diffs = zeros(length(test_ns));
fft_max_diffs = zeros(length(test_ns)); naive_max_diffs = zeros(length(test_ns));


for (i,n) in enumerate(test_ns)
    (theory_sq_output, naive_sq_output, fft_sq_output) = fft_n_points_test(n)
    fft_Δ = abs.(theory_sq_output .- fft_sq_output)
    fft_diffs[i] = sum(fft_Δ) / n
    fft_max_diffs[i] = maximum(fft_Δ)

    naive_Δ = abs.(theory_sq_output .- naive_sq_output)
    naive_diffs[i] = sum(naive_Δ) / n
    naive_max_diffs[i] = maximum(naive_Δ)
end

scatter(test_ns, [fft_diffs, naive_diffs], labels=[:fft, :naive], title="mean diffs",
    xlabel="# bins", ylabel="sum(abs(theory - numeric)) / # bins")

scatter(test_ns, [fft_max_diffs, naive_max_diffs], labels=[:fft, :naive], title="max diffs",
    xlabel="# bins", ylabel="max(abs(theory - numeric))")



test_ns = 71:10:1001
fft_diffs = zeros(length(test_ns)); naive_diffs = zeros(length(test_ns));
fft_max_diffs = zeros(length(test_ns)); naive_max_diffs = zeros(length(test_ns));


for (i,n) in enumerate(test_ns)
    (theory_sq_output, naive_sq_output, fft_sq_output) = fft_n_points_test(n)
    fft_Δ = abs.(theory_sq_output .- fft_sq_output)
    fft_diffs[i] = sum(fft_Δ) / n
    fft_max_diffs[i] = maximum(fft_Δ)

    naive_Δ = abs.(theory_sq_output .- naive_sq_output)
    naive_diffs[i] = sum(naive_Δ) / n
    naive_max_diffs[i] = maximum(naive_Δ)
end

scatter(test_ns, [fft_diffs, naive_diffs], labels=[:fft, :naive], title="mean diffs (odd)",
    xlabel="# bins", ylabel="sum(abs(theory - numeric)) / # bins")

scatter(test_ns, [fft_max_diffs, naive_max_diffs], labels=[:fft, :naive], title="max diffs (odd)",
    xlabel="# bins", ylabel="max(abs(theory - numeric))")


test_ns = 70:5:500
fft_diffs = zeros(length(test_ns)); naive_diffs = zeros(length(test_ns));
fft_max_diffs = zeros(length(test_ns)); naive_max_diffs = zeros(length(test_ns));


for (i,n) in enumerate(test_ns)
    (theory_sq_output, naive_sq_output, fft_sq_output) = fft_n_points_test(n)
    fft_Δ = abs.(theory_sq_output .- fft_sq_output)
    fft_diffs[i] = sum(fft_Δ) / n
    fft_max_diffs[i] = maximum(fft_Δ)

    naive_Δ = abs.(theory_sq_output .- naive_sq_output)
    naive_diffs[i] = sum(naive_Δ) / n
    naive_max_diffs[i] = maximum(naive_Δ)
end

scatter(test_ns, [fft_diffs, naive_diffs], labels=[:fft, :naive], title="mean diffs (5-step)",
    xlabel="# bins", ylabel="sum(abs(theory - numeric)) / # bins")

scatter(test_ns, [fft_max_diffs, naive_max_diffs], labels=[:fft, :naive], title="max diffs (5-step)",
    xlabel="# bins", ylabel="max(abs(theory - numeric))")
