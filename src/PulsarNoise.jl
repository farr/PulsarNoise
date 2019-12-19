module PulsarNoise

using DelimitedFiles

export read_tempo2_residuals

"""    read_tempo2_residuals(file)

Returns `(ts, rs, drs)`: times, residuals, and uncertainties that result from
reading `file`, a TEMPO2 residuals file.  In the event that there are multiple
residuals reported at a common time, the output will be a single residual that
is the appropriate weighted-average of the multiple values reported by TEMPO2.
"""
function read_tempo2_residuals(file)
    d = readdlm(file)

    inds = sortperm(d[:,1])
    ts = Float64[]
    rs = Float64[]
    drs = Float64[]

    tlast = d[inds[1],1]
    rlast = d[inds[1],2]
    drlast = d[inds[1],3]
    for (t, r, dr) in zip(d[inds[2:end],1], d[inds[2:end],2], d[inds[2:end],3])
        if t > tlast
            push!(ts, tlast)
            push!(rs, rlast)
            push!(drs, drlast)

            tlast = t
            rlast = r
            drlast = dr
        else
            wlast = 1.0/(drlast*drlast)
            wt = 1.0/(dr*dr)
            rlast = (wlast*rlast + wt*r)/(wlast+wt)
            drlast = sqrt(1.0/(wlast + wt))
        end
    end

    ts, rs, drs
end
end
