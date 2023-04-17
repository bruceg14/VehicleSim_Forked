# x is [p1 p2 v theata]
# u is [a angular velocity]
# Δ time step

"""
Unicycle model
"""
function f(position, quaternion, velocity, angular_vel, Δt)
    # r = angular_vel
    # mag = norm(r)

    # sᵣ = cos(mag*Δt / 2.0)
    # vᵣ = sin(mag*Δt / 2.0) * r/mag

    # sₙ = quaternion[1]
    # vₙ = quaternion[2:4]

    # s = sₙ*sᵣ - vₙ'*vᵣ
    # v = sₙ*vᵣ+sᵣ*vₙ+vₙ×vᵣ

    # new_position = position + Δt * velocity
    # new_quaternion = [s; v]
    # new_velocity = velocity
    # new_angular_vel = angular_vel
    # return [new_position; new_quaternion; new_velocity; new_angular_vel]

    r = angular_vel
    mag = norm(r)

    if mag < 1e-5
        sᵣ = 1.0
        vᵣ = zeros(3)
    else
        sᵣ = cos(mag*Δt / 2.0)
        vᵣ = sin(mag*Δt / 2.0) * axis
    end

    sₙ = quaternion[1]
    vₙ = quaternion[2:4]

    s = sₙ*sᵣ - vₙ'*vᵣ
    v = sₙ*vᵣ+sᵣ*vₙ+vₙ×vᵣ

    new_position = position + Δt * velocity
    new_quaternion = [s; v]
    new_velocity = velocity
    new_angular_vel = angular_vel
    return [new_position; new_quaternion; new_velocity; new_angular_vel]
end


"""
Jacobian of f_localization with respect to x, evaluated at x,Δ.
"""
function jac_fx(x, Δt)
    position = x[1:3]
    quaternion = x[4:7]
    velocity = x[8:10]
    angular_vel = x[11:13]
    r = angular_vel
    mag = norm(r)

    sᵣ = cos(mag*Δt / 2.0)
    vᵣ = sin(mag*Δt / 2.0) * r/mag

    sₙ = quaternion[1]
    vₙ = quaternion[2:4]

    s = sₙ*sᵣ - vₙ'*vᵣ
    v = sₙ*vᵣ+sᵣ*vₙ+vₙ×vᵣ

    # Some results I'm not sure about.
    # For Column 11
    jf11a = (4*(r[1]*conj(r[1]))^(1/2)*mag^2)
    jf11b = (2*(r[1]*conj(r[1]))^(1/2)*mag^2^(3/2))
    jf11c = abs(r[1])*(r[1] + conj(r[1]))

    jf_4_11 = (vᵣ[1]*mag*abs(r[1])*conj(vₙ[1])*(r[1] + conj(r[1])))/(2*(r[1]*conj(r[1]))^(1/2)*mag^2^(3/2)) - (sₙ*Δt*sin((Δt*mag)/2)*abs(r[1])*(r[1] + conj(r[1])))/(4*(r[1]*conj(r[1]))^(1/2)*mag) - (sin((Δt*mag)/2)*conj(vₙ[1]))/mag + (vᵣ[2]*mag*abs(r[1])*conj(vₙ[2])*(r[1] + conj(r[1])))/(2*(r[1]*conj(r[1]))^(1/2)*mag^2^(3/2)) + (vᵣ[3]*mag*abs(r[1])*conj(vₙ[3])*(r[1] + conj(r[1])))/(2*(r[1]*conj(r[1]))^(1/2)*mag^2^(3/2)) - (r[1]*Δt*sᵣ*abs(r[1])*conj(vₙ[1])*(r[1] + conj(r[1])))/(4*(r[1]*conj(r[1]))^(1/2)*mag^2) - (r[2]*Δt*sᵣ*abs(r[1])*conj(vₙ[2])*(r[1] + conj(r[1])))/(4*(r[1]*conj(r[1]))^(1/2)*mag^2) - (r[3]*Δt*sᵣ*abs(r[1])*conj(vₙ[3])*(r[1] + conj(r[1])))/(4*(r[1]*conj(r[1]))^(1/2)*mag^2)
    jf_5_11 = (sₙ*sin((Δt*mag)/2))/mag - (vₙ[1]*Δt*sin((Δt*mag)/2)*jf11c)/(4*(r[1]*conj(r[1]))^(1/2)*mag) - (r[1]*sₙ*sin((Δt*mag)/2)*jf11c)/jf11b + (r[2]*vₙ[3]*sin((Δt*mag)/2)*jf11c)/jf11b - (r[3]*vₙ[2]*sin((Δt*mag)/2)*jf11c)/jf11b + (r[1]*sₙ*Δt*sᵣ*jf11c)/jf11a - (r[2]*vₙ[3]*Δt*sᵣ*jf11c)/jf11a + (r[3]*vₙ[2]*Δt*sᵣ*jf11c)/jf11a
    jf_6_11 = (vₙ[3]*sin((Δt*mag)/2))/mag - (vₙ[2]*Δt*sin((Δt*mag)/2)*jf11c)/(4*(r[1]*conj(r[1]))^(1/2)*mag) - (r[2]*sₙ*sin((Δt*mag)/2)*jf11c)/jf11b - (r[1]*vₙ[3]*sin((Δt*mag)/2)*jf11c)/jf11b + (r[3]*vₙ[1]*sin((Δt*mag)/2)*jf11c)/jf11b + (r[2]*sₙ*Δt*sᵣ*jf11c)/jf11a + (r[1]*vₙ[3]*Δt*sᵣ*jf11c)/jf11a - (r[3]*vₙ[1]*Δt*sᵣ*jf11c)/jf11a
    jf_7_11 = (r[1]*vₙ[2]*sin((Δt*mag)/2)*jf11c)/jf11b - (vₙ[3]*Δt*sin((Δt*mag)/2)*jf11c)/(4*(r[1]*conj(r[1]))^(1/2)*mag) - (vₙ[2]*sin((Δt*mag)/2))/mag - (r[2]*vₙ[1]*sin((Δt*mag)/2)*jf11c)/jf11b - (r[3]*sₙ*sin((Δt*mag)/2)*jf11c)/jf11b - (r[1]*vₙ[2]*Δt*sᵣ*jf11c)/jf11a + (r[2]*vₙ[1]*Δt*sᵣ*jf11c)/jf11a + (r[3]*sₙ*Δt*sᵣ*jf11c)/jf11a

    # For Column 12
    jf12a = (4*(r[2]*conj(r[2]))^(1/2)*mag^2)
    jf12b = (2*(r[2]*conj(r[2]))^(1/2)*mag^2^(3/2))
    jf12c = abs(r[2])*(r[2] + conj(r[2]))

    jf_4_12 = (vᵣ[1]*mag*abs(r[2])*conj(vₙ[1])*(r[2] + conj(r[2])))/jf12b - (sₙ*Δt*sin((Δt*mag)/2)*jf12c)/(4*(r[2]*conj(r[2]))^(1/2)*mag) - (sin((Δt*mag)/2)*conj(vₙ[2]))/mag + (vᵣ[2]*mag*abs(r[2])*conj(vₙ[2])*(r[2] + conj(r[2])))/jf12b + (vᵣ[3]*mag*abs(r[2])*conj(vₙ[3])*(r[2] + conj(r[2])))/jf12b - (r[1]*Δt*sᵣ*abs(r[2])*conj(vₙ[1])*(r[2] + conj(r[2])))/jf12a - (r[2]*Δt*sᵣ*abs(r[2])*conj(vₙ[2])*(r[2] + conj(r[2])))/jf12a - (r[3]*Δt*sᵣ*abs(r[2])*conj(vₙ[3])*(r[2] + conj(r[2])))/jf12a
    jf_5_12 = (r[2]*vₙ[3]*sin((Δt*mag)/2)*jf12c)/jf12b - (vₙ[1]*Δt*sin((Δt*mag)/2)*jf12c)/(4*(r[2]*conj(r[2]))^(1/2)*mag) - (r[1]*sₙ*sin((Δt*mag)/2)*jf12c)/jf12b - (vₙ[3]*sin((Δt*mag)/2))/mag - (r[3]*vₙ[2]*sin((Δt*mag)/2)*jf12c)/jf12b + (r[1]*sₙ*Δt*sᵣ*jf12c)/jf12a - (r[2]*vₙ[3]*Δt*sᵣ*jf12c)/jf12a + (r[3]*vₙ[2]*Δt*sᵣ*jf12c)/jf12a
    jf_6_12 = (sₙ*sin((Δt*mag)/2))/mag - (vₙ[2]*Δt*sin((Δt*mag)/2)*jf12c)/(4*(r[2]*conj(r[2]))^(1/2)*mag) - (r[2]*sₙ*sin((Δt*mag)/2)*jf12c)/jf12b - (r[1]*vₙ[3]*sin((Δt*mag)/2)*jf12c)/jf12b + (r[3]*vₙ[1]*sin((Δt*mag)/2)*jf12c)/jf12b + (r[2]*sₙ*Δt*sᵣ*jf12c)/jf12a + (r[1]*vₙ[3]*Δt*sᵣ*jf12c)/jf12a - (r[3]*vₙ[1]*Δt*sᵣ*jf12c)/jf12a
    jf_7_12 = (vₙ[1]*sin((Δt*mag)/2))/mag - (vₙ[3]*Δt*sin((Δt*mag)/2)*jf12c)/(4*(r[2]*conj(r[2]))^(1/2)*mag) + (r[1]*vₙ[2]*sin((Δt*mag)/2)*jf12c)/jf12b - (r[2]*vₙ[1]*sin((Δt*mag)/2)*jf12c)/jf12b - (r[3]*sₙ*sin((Δt*mag)/2)*jf12c)/jf12b - (r[1]*vₙ[2]*Δt*sᵣ*jf12c)/jf12a + (r[2]*vₙ[1]*Δt*sᵣ*jf12c)/jf12a + (r[3]*sₙ*Δt*sᵣ*jf12c)/jf12a

    # For Column 13
    jf13a = (4*(r[3]*conj(r[3]))^(1/2)*mag^2)
    jf13b = (2*(r[3]*conj(r[3]))^(1/2)*mag^2^(3/2))
    jf13c = abs(r[3])*(r[3] + conj(r[3]))

    jf_4_13 = (vᵣ[1]*mag*abs(r[3])*conj(vₙ[1])*(r[3] + conj(r[3])))/jf13b - (sₙ*Δt*sin((Δt*mag)/2)*jf13c)/(4*(r[3]*conj(r[3]))^(1/2)*mag) - (sin((Δt*mag)/2)*conj(vₙ[3]))/mag + (vᵣ[2]*mag*abs(r[3])*conj(vₙ[2])*(r[3] + conj(r[3])))/jf13b + (vᵣ[3]*mag*abs(r[3])*conj(vₙ[3])*(r[3] + conj(r[3])))/jf13b - (r[1]*Δt*sᵣ*abs(r[3])*conj(vₙ[1])*(r[3] + conj(r[3])))/jf13a - (r[2]*Δt*sᵣ*abs(r[3])*conj(vₙ[2])*(r[3] + conj(r[3])))/jf13a - (r[3]*Δt*sᵣ*abs(r[3])*conj(vₙ[3])*(r[3] + conj(r[3])))/jf13a
    jf_5_13 = (vₙ[2]*sin((Δt*mag)/2))/mag - (vₙ[1]*Δt*sin((Δt*mag)/2)*jf13c)/(4*(r[3]*conj(r[3]))^(1/2)*mag) - (r[1]*sₙ*sin((Δt*mag)/2)*jf13c)/jf13b + (r[2]*vₙ[3]*sin((Δt*mag)/2)*jf13c)/jf13b - (r[3]*vₙ[2]*sin((Δt*mag)/2)*jf13c)/jf13b + (r[1]*sₙ*Δt*sᵣ*jf13c)/jf13a - (r[2]*vₙ[3]*Δt*sᵣ*jf13c)/jf13a + (r[3]*vₙ[2]*Δt*sᵣ*jf13c)/jf13a
    jf_6_13 = (r[3]*vₙ[1]*sin((Δt*mag)/2)*jf13c)/jf13b - (vₙ[2]*Δt*sin((Δt*mag)/2)*jf13c)/(4*(r[3]*conj(r[3]))^(1/2)*mag) - (r[2]*sₙ*sin((Δt*mag)/2)*jf13c)/jf13b - (r[1]*vₙ[3]*sin((Δt*mag)/2)*jf13c)/jf13b - (vₙ[1]*sin((Δt*mag)/2))/mag + (r[2]*sₙ*Δt*sᵣ*jf13c)/jf13a + (r[1]*vₙ[3]*Δt*sᵣ*jf13c)/jf13a - (r[3]*vₙ[1]*Δt*sᵣ*jf13c)/jf13a
    jf_7_13 = (sₙ*sin((Δt*mag)/2))/mag - (vₙ[3]*Δt*sin((Δt*mag)/2)*jf13c)/(4*(r[3]*conj(r[3]))^(1/2)*mag) + (r[1]*vₙ[2]*sin((Δt*mag)/2)*jf13c)/jf13b - (r[2]*vₙ[1]*sin((Δt*mag)/2)*jf13c)/jf13b - (r[3]*sₙ*sin((Δt*mag)/2)*jf13c)/jf13b - (r[1]*vₙ[2]*Δt*sᵣ*jf13c)/jf13a + (r[2]*vₙ[1]*Δt*sᵣ*jf13c)/jf13a + (r[3]*sₙ*Δt*sᵣ*jf13c)/jf13a


    [1. 0. 0.   0.      0.      0.      0.      Δt 0. 0. 0.         0.          0.;
     0. 1. 0.   0.      0.      0.      0.      0. Δt 0. 0.         0.          0.;
     0. 0. 1.   0.      0.      0.      0.      0. 0. Δt 0.         0.          0.;
     0. 0. 0.   sᵣ      -vᵣ[1]  -vᵣ[2]  -vᵣ[3]  0. 0. 0. jf_4_11    jf_4_12     jf_4_13;
     0. 0. 0.   vᵣ[1]   sᵣ      vᵣ[3]   -vᵣ[2]  0. 0. 0. jf_5_11    jf_5_12     jf_5_13;
     0. 0. 0.   vᵣ[2]   -vᵣ[3]  sᵣ      vᵣ[1]   0. 0. 0. jf_6_11    jf_6_12     jf_6_13;
     0. 0. 0.   vᵣ[3]   vᵣ[2]   -vᵣ[1]  sᵣ      0. 0. 0. jf_7_11    jf_7_12     jf_7_13;
     0. 0. 0.   0.      0.      0.      0.      1. 0. 0. 0.         0.          0.;
     0. 0. 0.   0.      0.      0.      0.      0. 1. 0. 0.         0.          0.;
     0. 0. 0.   0.      0.      0.      0.      0. 0. 1. 0.         0.          0.;
     0. 0. 0.   0.      0.      0.      0.      0. 0. 0. 1.         0.          0.;
     0. 0. 0.   0.      0.      0.      0.      0. 0. 0. 0.         1.          0.;
     0. 0. 0.   0.      0.      0.      0.      0. 0. 0. 0.         0.          1.]

end



"""
Jacobian of f with respect to u, evaluated at x,u,Δ.
"""
# function jac_fu(x, u, Δ)
#     [0 0;
#      0 0;
#      Δ 0;
#      0 Δ]
# end

"""
Non-standard measurement model. Can we extract state estimate from just this?
"""
# x is position, quaternion, velocity, angular vel
function h_gps(x)
    
    T = get_gps_transform()
    gps_loc_body = T*[zeros(3); 1.0]


    xyz_body = x[1:3]
    q_body = x[4:7]
    

    Tbody = get_body_transform(q_body, xyz_body)
    xyz_gps = Tbody * [gps_loc_body; 1]
    meas = xyz_gps[1:2] 
    gps_meas = GPSMeasurement(t, meas)

    return meas
end

function h_imu(x)
    velocity = x[8: 10]
    angular_vel = [11:13]
    # Do Imu now
    T_body_imu = get_imu_transform()
    T_imu_body = invert_transform(T_body_imu)
    R_imu = T_imu_body[1:3,1:3]
    p_imu = T_imu_body[1:3,end]

    w_imu = R_imu * angular_vel
    v_imu = R_imu * velocity + p_imu × w_imu


    return [v_imu, w_imu]
end

"""
Jacobian of h with respect to x, evaluated at x.
"""
function jac_h_gps(x)
    [1.0 0.0 0.0    (26*x[6])/5-6*x[4]-2*x[7]       2*x[6]-6*x[5]+(26*x[7])/5       (26*x[4])/5+2*x[5]+6*x[6]       (26*x[5])/5-2*x[4]+6*x[7]       0.0     0.0     0.0     0.0     0.0     0.0;
     0.0 1.0 0.0    2*x[4]-(26*x[5])/5-6*x[7]       -(26*x[4])/5-2*x[5]-6*x[6]      2*x[6]-6*x[5]+(26*x[7])/5       (26*x[6])/5-6*x[4]-2*x[7]       0.0     0.0     0.0     0.0     0.0     0.0;]
end


    

function jac_h_imu(x)
   [0.0 0.0 0.0    0.0                             0.0                             0.0                             0.0                             0.9998  0.0     -0.0199 0.0     0.7     0.0;
    0.0 0.0 0.0    0.0                             0.0                             0.0                             0.0                             0.0     1.0     0.0     -0.7    0.0     0.0;
    0.0 0.0 0.0    0.0                             0.0                             0.0                             0.0                             0.0199  0.0     0.9998  0.0     0.014   0.0;
    0.0 0.0 0.0    0.0                             0.0                             0.0                             0.0                             0.0     0.0     0.0     0.9998  0.0     -0.0199;
    0.0 0.0 0.0    0.0                             0.0                             0.0                             0.0                             0.0     0.0     0.0     0.0     1.0     0.0;
    0.0 0.0 0.0    0.0                             0.0                             0.0                             0.0                             0.0     0.0     0.0     0.0199  0.0     0.9998]
end



"""
Extended kalman filter implementation.

Assume that the 'true' physical update in the world is given by 

xₖ = f(xₖ₋₁, uₖ, Δ), where Δ is the time difference between times k and k-1.

Here, uₖ is the 'true' controls applied to the system. These controls can be assumed to be a random variable,
with probability distribution given by 𝒩 (mₖ, proc_cov) where mₖ is some IMU-like measurement, and proc_cov is a constant covariance matrix.

The process model distribution is then approximated as:

P(xₖ | xₖ₋₁, uₖ) ≈ 𝒩 ( Axₖ₋₁  + c, Σ̂ )

where 
A = ∇ₓf(μₖ₋₁,  Δ),

c = f(μₖ₋₁, Δ) - Aμₖ₋₁ 

μ̂ = Aμₖ₋₁  + c
  = f(μₖ₋₁, mₖ, Δ)
Σ̂ = A Σₖ₋₁ A' + proc_cov(how much noisy is for the process model), 


Further, assume that the 'true' measurement generation in the world is given by

zₖ = h(xₖ) + wₖ,

where wₖ is some additive gaussian noise with probability density function given by

𝒩 (0, meas_var).

The measurement model is then approximated as 

P(zₖ | xₖ) ≈ 𝒩 ( C xₖ + d , meas_var )


where 
C = ∇ₓ h(μ̂), 
d = h(μ̂) - Cμ̂

The extended Kalman filter update equations can be implemented as the following:

Σₖ = (Σ̂⁻¹ + C' (meas_var)⁻¹ C)⁻¹
μₖ = Σₖ ( Σ̂⁻¹ μ̂ + C' (meas_var)⁻¹ (zₖ - d) )

"""
function filter(meas; Δ,  μK , Σk, localization_state_channel)
    meas_var = Diagonal([1,1,0.05,0.1,0.1,0.1,0.1,0.05,0.05,0.05,0.05,0.05,0.05]) 
    proc_cov = Diagonal([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

    x_prev = μK[μK.length]

    Σ_prev = Σk[Σk.length]

    xₖ = f(x_prev[1:3], x_prev[4:7], x_prev[8:10], x_prev[11:13], Δ) #Current x


    z = meas

    A = jac_fx(xₖ, Δ)
    
    μ̂ = f(x_prev,  Δ)

    Σ̂ = A*Σ_prev*A' + proc_cov

    if meas isa GPSMeasurement
        C = jac_h_gps(μ̂) # gps version
        d = h_gps(μ̂) - C*μ̂ # gps version
        Σ = inv(inv(Σ̂) + C'*inv(meas_var)*C)
        μ = Σ * (inv(Σ̂) * μ̂ + C'*inv(meas_var) * (z - d))

        full_state = FullVehicleState{position : μ[1:3], velocity: μ[8:10], ang_vel:[11:13], orientation:[4:7]}

        local_state = MyLocalizationType{last_update: meas.time, x : full_state}
        put!(localization_state_channel, local_state)
    else
        C = jac_h_imu(μ̂) # imu version
        d = h_imu(μ̂) - C*μ̂ # imu version
        Σ = inv(inv(Σ̂) + C'*inv(meas_var)*C)
        μ = Σ * (inv(Σ̂) * μ̂ + C'*inv(meas_var) * (z - d))
        
        full_state = FullVehicleState{position : μ[1:3], velocity: μ[8:10], ang_vel:[11:13], orientation:[4:7]}

        local_state = MyLocalizationType{last_update: meas.time, x : full_state}

        put!(localization_state_channel, local_state)
    end


end

struct FullVehicleState 
    position::SVector{3, Float64}
    velocity::SVector{3, Float64}
    ang_vel::SVector{3, Float64}
    orientation::SVector{3, Float64}
end

struct MyLocalizationType
    last_update:: Float64
    x::FullVehicleState
end