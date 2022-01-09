using DelimitedFiles
function loadcurve( filename )
    ## source from: http://www.pjm.com/markets-and-operations/ops-analysis/historical-load-data.aspx
    temp_data = readdlm( filename, ',',header=true)[1][:,2:end]
    raw_data = Array{Float64}(undef,size(temp_data))
    normalized_data = Array{Float64}(undef,size(temp_data))
    raw_data[:] = temp_data[:]
    N_days = size(raw_data,1)
    for k = 1:N_days
        normalized_data[k,:] = raw_data[k,:] / reduce(max,raw_data[k,:])
    end
    load24 = normalized_data
    return load24
end
