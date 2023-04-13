struct FullVehicleState 
    time::Float64
    position::SVector{3, Float64}
    velocity::SVector{3, Float64}
    ang_vel::SVector{3, Float64}
    orientation::SVector{3, Float64}
end

function localization(gps_channel, imu_channel, localization_state_channel)
    μ0 =[0,0, 2.645, zeros(10)]
    μk= [μ0,]

    Σ=Diagonal([50,50,0.5,1, 1, 1, 1, 3, 3, 3, 0.5, 0.5, 0.5])
    Σk = Matrix{Float64}[Σ,]
    

    initial_gps = {time: 0, lat: 0, long: 0}
    initial_imu = {time: 0, linear_vel : [0,0,0], angular_vel: [0,0,0]}
    fresh_gps_meas = [initial_gps,] # Just to set the initial time to be 0
    fresh_imu_meas = [initial_imu,]
    latest_meas_time = -Inf
    while true
        fresh_gps_meas = [] # Just to set the initial time to be 0
        fresh_imu_meas = []
        all_meas = []
        while isready(gps_channel)
            meas = take!(gps_channel)
            push!(all_meas, meas)
        end
        while isready(imu_channel)
            meas = take!(imu_channel)
            push!(all_meas, meas)
        end

        sorted_all_meas = sort(all_meas, by = x -> x.time)

        start = 1
        for(i = 1; i <= sorted_all_meas.length; ++i){
            if(sorted_all_meas[i].time >= latest_meas_time){
                start = i 
                break
            }
        }
        # throw away any measuremetns in all_meas that are from BEFORE last_meas_time

        for (i = start, i <=sorted_all_meas.length; ++i)
            dt = sorted_all_meas[i].time - latest_meas_time
            filter(sorted_all_meas[i], dt, μK , Σk,  localization_state_channel)
            latest_meas_time = meas.time
        end

    end
end



        



