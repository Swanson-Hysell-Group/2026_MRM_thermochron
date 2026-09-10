# Export the posterior time-temperature path density of a Thermochron.jl
# inversion as an 8-bit array, so that composite figures can be assembled
# outside Julia (see ../code/make_inversion_figures.py).
#
# The density image is computed exactly as in tTinversion_mrm.jl (1800 x 1200
# bins spanning the full model time and temperature range), normalised to its
# 98.5th percentile and clipped to [0, 1] -- the same scaling that
# tTinversion_mrm.jl applies before drawing `<prefix>_tT.pdf` -- and written
# as an unsigned 8-bit NumPy `.npy` array of shape (1200, 1800): rows run from
# 0 to Tinit (°C, top to bottom), columns from 0 to tinit (Ma, left to
# right). Time and temperature ranges are recorded in the companion
# `<prefix>_tT_density.json`.
#
# The full posteriors (`.jld2`, ~800 MB each) are not archived in the
# repository, so this script can only be run on the machine that holds them.
# The exported density maps are converted to PNG and archived by
# make_inversion_figures.py.
#
# Usage, from the inversion_files/ directory, with the batch's pinned
# environment:
#
#   julia --project=swansonhysell_runs_2026-08a export_tT_density.jl \
#       swansonhysell_runs_2026-08a/MRM_all_constraints_seed20260710 [...]
#
# Each argument is an output prefix (the .jld2 path without its extension).

using Thermochron
using JLD2, CodecZlib   # nanpctile and nanmean are re-exported by Thermochron
cd(@__DIR__)

"""
    write_npy_uint8(path, A)

Write a matrix of UInt8 as a NumPy `.npy` (format version 1.0) file in C
order, so that `numpy.load(path)[i, j] == A[i, j]`.
"""
function write_npy_uint8(path::AbstractString, A::AbstractMatrix{UInt8})
    header = "{'descr': '|u1', 'fortran_order': False, 'shape': ($(size(A,1)), $(size(A,2))), }"
    # magic (6) + version (2) + header length (2) + header + newline must be a
    # multiple of 64 bytes
    pad = (64 - (10 + length(header) + 1) % 64) % 64
    header *= " "^pad * "\n"
    open(path, "w") do io
        write(io, UInt8(0x93)); write(io, "NUMPY")
        write(io, UInt8(1)); write(io, UInt8(0))
        write(io, htol(UInt16(length(header))))
        write(io, header)
        write(io, permutedims(A))   # column-major of the transpose == C order of A
    end
end

for prefix in ARGS
    jld = prefix * ".jld2"
    isfile(jld) || error("Posterior file not found: $jld")
    @info "Loading $jld"
    @time @load jld tT model
    tinit, Tinit, Tnow, tnow = model.tinit, model.Tinit, model.Tnow, model.tnow

    # Same call as in tTinversion_mrm.jl
    @time (tTimage, xc, yc) = image_from_paths!(tT; xresolution=1800, yresolution=1200,
        xrange=[tnow, tinit], yrange=[Tnow, Tinit])
    @info "Density image size $(size(tTimage)); x bins $(first(xc)) to $(last(xc)) Ma, y bins $(first(yc)) to $(last(yc)) °C"

    # Normalise as in tTinversion_mrm.jl: imsc(tTimage, ylcn, 0, nanpctile(tTimage[:], 98.5))
    scale = nanpctile(vec(tTimage), 98.5)
    A = UInt8.(round.(clamp.(tTimage ./ scale, 0, 1) .* 255))
    @info "98.5th-percentile normalisation: $scale"

    out = prefix * "_tT_density.npy"
    write_npy_uint8(out, A)
    open(prefix * "_tT_density.json", "w") do io
        write(io, """{"source": "$(basename(jld))", "shape": [$(size(A,1)), $(size(A,2))],
 "time_Ma": [$(tnow), $(tinit)], "temperature_C": [$(Tnow), $(Tinit)],
 "normalisation": "98.5th percentile of the path-density image, clipped to [0, 1] and scaled to 0-255",
 "post_burnin_steps": $(size(tT.Tpointdist, 2)), "burnin": $(model.burnin),
 "mean_loglikelihood": $(nanmean(view(tT.lldist, model.burnin:model.nsteps))),
 "mean_acceptance": $(nanmean(view(tT.acceptancedist, model.burnin:model.nsteps))),
 "mean_npoints": $(nanmean(view(tT.ndist, model.burnin:model.nsteps)))}
""")
    end
    @info "Wrote $out"
    tT = nothing; GC.gc()
end
