from dataclasses import dataclass
from math import *
import matplotlib.pyplot as plt
import random
import statistics

@dataclass
class Rocket:
    dry_mass_kg: float
    diameter_m: float
    length_m: float
    mach_points: list[float]
    cd_points: list[float]
    cn_alpha: float
    cm_alpha: float
    moi_kg_m2: float

    @property
    def frontal_area_m2(self) -> float:
        radius = self.diameter_m / 2.0
        area = pi * radius**2
        return area

@dataclass
class Motor:
    thrust_N: float
    burn_time_s: float
    prop_mass_kg: float

    @property
    def mass_flow_rate_kg_s(self) -> float:
        if self.burn_time_s <= 0:
            raise ValueError("burn_time_s must be positive")

        mdot = self.prop_mass_kg / self.burn_time_s
        return mdot

@dataclass
class Vehicle:
    rocket: Rocket
    motor: Motor

    @property
    def initial_mass_kg(self) -> float:
        total_mass = self.rocket.dry_mass_kg + self.motor.prop_mass_kg
        return total_mass

@dataclass
class Launch:
    angle_deg: float

@dataclass
class Parachute:
    cd: float
    area_m2: float
    deploy_delay_s: float
    rot_damping_1_s: float
    rot_stiffness_1_s2: float
    theta_eq_rad: float

@dataclass
class Environment:
    wind_altitudes_m: list[float]
    wind_speeds_m_s: list[float]

###########################################################################################
def scaled_environment(environment: Environment, scale: float) -> Environment:
    wind_altitudes_m = environment.wind_altitudes_m[:]

    wind_speeds_m_s = []
    for i in range(len(environment.wind_speeds_m_s)):
        val = environment.wind_speeds_m_s[i] * scale
        wind_speeds_m_s.append(val)

    env_scaled = Environment(
        wind_altitudes_m=wind_altitudes_m,
        wind_speeds_m_s=wind_speeds_m_s
    )

    return env_scaled

###########################################################################################
def cd_from_mach(mach: float, rocket: Rocket) -> float:
    ### Inputs ###
    mach_pts = rocket.mach_points
    cd_pts = rocket.cd_points

    ### Checks ###
    if len(mach_pts) != len(cd_pts):
        raise ValueError("mach_points and cd_points must have the same length")

    if len(mach_pts) < 2:
        raise ValueError("Need at least two Mach/CD points")

    for i in range(len(mach_pts) - 1):
        if mach_pts[i] >= mach_pts[i + 1]:
            raise ValueError("mach_points must be strictly increasing")

    ### Endpoints ###
    if mach <= mach_pts[0]:
        return cd_pts[0]

    if mach >= mach_pts[-1]:
        return cd_pts[-1]

    ### Linear interpolation ###
    for i in range(len(mach_pts) - 1):
        m0 = mach_pts[i]
        m1 = mach_pts[i + 1]
        cd0 = cd_pts[i]
        cd1 = cd_pts[i + 1]

        if m0 <= mach <= m1:
            frac = (mach - m0) / (m1 - m0)
            cd_val = cd0 + frac * (cd1 - cd0)
            return cd_val

    return cd_pts[-1]

###########################################################################################
def thrust_force(t: float, motor: Motor) -> float:
    if 0.0 <= t <= motor.burn_time_s:
        return motor.thrust_N

    return 0.0

###########################################################################################
def gravity_accel(alt_m: float) -> float:
    g0 = 9.80665
    r_earth = 6371000.0

    ratio = r_earth / (r_earth + alt_m)
    g = g0 * ratio**2

    return g

###########################################################################################
def signed_gravity_force(mass_kg: float, alt_m: float) -> float:
    g = gravity_accel(alt_m)
    gravity_force = -mass_kg * g
    return gravity_force

###########################################################################################
def drag_force_magnitude(z_m: float, vx_m_s: float, vz_m_s: float, rocket: Rocket,
    environment: Environment, atmosphere_func) -> tuple[float, float, float, float]:
    ### Altitude ###
    z_m = max(0.0, z_m)

    ### Relative velocity ###
    wind_x = environment.wind_speed_m_s
    vrel_x = vx_m_s - wind_x
    vrel_z = vz_m_s

    vrel_sq = vrel_x**2 + vrel_z**2
    v_rel = sqrt(vrel_sq)

    ### Atmosphere ###
    _, _, rho, a = atmosphere_si(z_m, atmosphere_func)

    ### Aero values ###
    if a > 0:
        mach = v_rel / a
    else:
        mach = 0.0

    cd = cd_from_mach(mach, rocket)
    q = 0.5 * rho * v_rel**2
    drag_mag = q * cd * rocket.frontal_area_m2

    return drag_mag, mach, cd, q

###########################################################################################
def drag_force_components(z_m: float, vx_m_s: float, vz_m_s: float, rocket: Rocket,
    environment: Environment, atmosphere_func) -> tuple[float, float, float, float, float]:
    ### Altitude ###
    z_m = max(0.0, z_m)

    ### Relative velocity ###
    wind_x = environment.wind_speed_m_s
    vrel_x = vx_m_s - wind_x
    vrel_z = vz_m_s

    vrel_sq = vrel_x**2 + vrel_z**2
    v_rel = sqrt(vrel_sq)

    if v_rel == 0.0:
        cd = cd_from_mach(0.0, rocket)
        return 0.0, 0.0, 0.0, cd, 0.0

    ### Drag magnitude ###
    drag_mag, mach, cd, q = drag_force_magnitude(
        z_m, vx_m_s, vz_m_s, rocket, environment, atmosphere_func
    )

    drag_x = -drag_mag * (vrel_x / v_rel)
    drag_z = -drag_mag * (vrel_z / v_rel)

    return drag_x, drag_z, mach, cd, q

###########################################################################################
def aero_forces_and_moment(z_m: float, vx_m_s: float, vz_m_s: float,
    theta_rad: float, rocket: Rocket, parachute: Parachute,
    environment: Environment, atmosphere_func, parachute_deployed: bool):
    ### Altitude ###
    z_m = max(0.0, z_m)

    ### Relative velocity ###
    wind_x = wind_velocity_x(z_m, environment)
    vrel_x = vx_m_s - wind_x
    vrel_z = vz_m_s

    vrel_sq = vrel_x**2 + vrel_z**2
    v_rel = sqrt(vrel_sq)

    ### Atmosphere ###
    _, _, rho, a = atmosphere_si(z_m, atmosphere_func)

    if a > 0:
        mach = v_rel / a
    else:
        mach = 0.0

    q = 0.5 * rho * v_rel**2

    ### Relative flow direction ###
    if v_rel > 0.0:
        ux = vrel_x / v_rel
        uz = vrel_z / v_rel
        gamma = atan2(vrel_z, vrel_x)
    else:
        ux = 1.0
        uz = 0.0
        gamma = theta_rad

    ### Parachute descent ###
    if parachute_deployed:
        cd = parachute.cd
        area = parachute.area_m2
        drag_mag = q * cd * area

        drag_x = -drag_mag * ux
        drag_z = -drag_mag * uz

        aero_x = drag_x
        aero_z = drag_z
        moment = 0.0
        alpha = 0.0

        return aero_x, aero_z, moment, mach, cd, q, alpha, v_rel

    ### Rocket body aerodynamics ###
    cd = cd_from_mach(mach, rocket)
    drag_mag = q * cd * rocket.frontal_area_m2

    alpha = theta_rad - gamma

    cn = rocket.cn_alpha * alpha
    normal_mag = q * cn * rocket.frontal_area_m2

    cm = rocket.cm_alpha * alpha
    moment = q * cm * rocket.frontal_area_m2 * rocket.length_m

    drag_x = -drag_mag * ux
    drag_z = -drag_mag * uz

    normal_x = -normal_mag * uz
    normal_z = normal_mag * ux

    aero_x = drag_x + normal_x
    aero_z = drag_z + normal_z

    return aero_x, aero_z, moment, mach, cd, q, alpha, v_rel

###########################################################################################
def signed_thrust_force(t: float, motor: Motor) -> float:
    thrust = thrust_force(t, motor)
    return thrust

###########################################################################################
def mass_rate(t: float, vehicle: Vehicle, mass_kg: float) -> float:
    dry_mass = vehicle.rocket.dry_mass_kg

    if 0.0 <= t <= vehicle.motor.burn_time_s and mass_kg > dry_mass:
        mdot = -vehicle.motor.mass_flow_rate_kg_s
        return mdot

    return 0.0

###########################################################################################
def state_derivatives(t: float, state: list[float], vehicle: Vehicle,
    launch: Launch, parachute: Parachute, environment: Environment,
    atmosphere_func, parachute_deployed: bool):
    ### State unpacking ###
    x_m = state[0]
    z_m = state[1]
    vx_m_s = state[2]
    vz_m_s = state[3]
    theta_rad = state[4]
    omega_rad_s = state[5]
    mass_kg = state[6]

    z_m = max(0.0, z_m)

    ### Aerodynamic forces and moments ###
    aero_x, aero_z, moment, _, _, _, _, _ = aero_forces_and_moment(
        z_m,
        vx_m_s,
        vz_m_s,
        theta_rad,
        vehicle.rocket,
        parachute,
        environment,
        atmosphere_func,
        parachute_deployed,
    )

    ### Gravity ###
    gravity_z = signed_gravity_force(mass_kg, z_m)

    ### Thrust, mass change, and rotational dynamics ###
    if parachute_deployed:
        thrust_x = 0.0
        thrust_z = 0.0
        dmdt = 0.0

        theta_diff = theta_rad - parachute.theta_eq_rad
        sin_term = sin(theta_diff)
        cos_term = cos(theta_diff)
        angle_error = atan2(sin_term, cos_term)

        dthetadt = omega_rad_s
        domegadt = -parachute.rot_stiffness_1_s2 * angle_error \
            - parachute.rot_damping_1_s * omega_rad_s
    else:
        thrust_mag = signed_thrust_force(t, vehicle.motor)
        thrust_x = thrust_mag * cos(theta_rad)
        thrust_z = thrust_mag * sin(theta_rad)
        dmdt = mass_rate(t, vehicle, mass_kg)

        dthetadt = omega_rad_s
        domegadt = moment / vehicle.rocket.moi_kg_m2

    ### Translational dynamics ###
    force_x = thrust_x + aero_x
    force_z = thrust_z + aero_z + gravity_z

    dxdt = vx_m_s
    dzdt = vz_m_s
    dvxdt = force_x / mass_kg
    dvzdt = force_z / mass_kg

    derivs = [dxdt, dzdt, dvxdt, dvzdt, dthetadt, domegadt, dmdt]
    return derivs

###########################################################################################
def rk4_step(f, t, state, dt, vehicle, launch, parachute,
    environment, atmosphere_func, parachute_deployed):
    ### k1 ###
    k1 = f(t, state, vehicle, launch, parachute, environment,
        atmosphere_func, parachute_deployed)

    ### k2 ###
    k2_state = []
    for i in range(len(state)):
        val = state[i] + 0.5 * dt * k1[i]
        k2_state.append(val)

    k2 = f(t + 0.5 * dt, k2_state, vehicle, launch, parachute,
        environment, atmosphere_func, parachute_deployed)

    ### k3 ###
    k3_state = []
    for i in range(len(state)):
        val = state[i] + 0.5 * dt * k2[i]
        k3_state.append(val)

    k3 = f(t + 0.5 * dt, k3_state, vehicle, launch, parachute,
        environment, atmosphere_func, parachute_deployed)

    ### k4 ###
    k4_state = []
    for i in range(len(state)):
        val = state[i] + dt * k3[i]
        k4_state.append(val)

    k4 = f(t + dt, k4_state, vehicle, launch, parachute,
        environment, atmosphere_func, parachute_deployed)

    ### Final RK4 update ###
    new_state = []
    for i in range(len(state)):
        incr = (dt / 6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])
        val = state[i] + incr
        new_state.append(val)

    return new_state

###########################################################################################
def atmosphere_si(alt_m: float, atmosphere_func):
    alt_ft = alt_m * 3.28084
    atm_vals = atmosphere_func(alt_ft)
    return atm_vals

###########################################################################################
def mach_number(alt_m: float, vel_m_s: float, atmosphere_func) -> float:
    _, _, _, a = atmosphere_si(alt_m, atmosphere_func)

    if a > 0:
        mach = abs(vel_m_s) / a
    else:
        mach = 0.0

    return mach

###########################################################################################
def dynamic_pressure(alt_m: float, vel_m_s: float, atmosphere_func) -> float:
    _, _, rho, _ = atmosphere_si(alt_m, atmosphere_func)
    q = 0.5 * rho * vel_m_s**2
    return q

###########################################################################################
def wind_velocity_x(alt_m: float, environment: Environment) -> float:
    """
    Returns horizontal wind speed [m/s] interpolated from altitude table.
    """
    ### Profile data ###
    alts = environment.wind_altitudes_m
    winds = environment.wind_speeds_m_s

    ### Checks ###
    if len(alts) != len(winds):
        raise ValueError("wind_altitudes_m and wind_speeds_m_s must have the same length")

    if len(alts) < 2:
        raise ValueError("Need at least two wind profile points")

    for i in range(len(alts) - 1):
        if alts[i] >= alts[i + 1]:
            raise ValueError("wind_altitudes_m must be strictly increasing")

    ### Endpoints ###
    if alt_m <= alts[0]:
        return winds[0]

    if alt_m >= alts[-1]:
        return winds[-1]

    ### Linear interpolation ###
    for i in range(len(alts) - 1):
        z0 = alts[i]
        z1 = alts[i + 1]
        w0 = winds[i]
        w1 = winds[i + 1]

        if z0 <= alt_m <= z1:
            frac = (alt_m - z0) / (z1 - z0)
            wind_val = w0 + frac * (w1 - w0)
            return wind_val

    return winds[-1]

###########################################################################################
def simulate(vehicle: Vehicle, launch: Launch, parachute: Parachute,
    environment: Environment, atmosphere_func, dt=0.01, t_max=300.0):
    ### Initial conditions ###
    t = 0.0
    x = 0.0
    z = 0.0
    vx = 0.0
    vz = 0.0
    theta = radians(launch.angle_deg)
    omega = 0.0
    mass = vehicle.initial_mass_kg

    state = [x, z, vx, vz, theta, omega, mass]

    ### Histories ###
    time_hist = []
    x_hist = []
    z_hist = []
    vx_hist = []
    vz_hist = []
    speed_hist = []
    theta_hist = []
    omega_hist = []
    alpha_hist = []
    mass_hist = []
    accel_hist = []
    mach_hist = []
    q_hist = []
    cd_hist = []

    ### Events / metrics ###
    burnout_time = None
    burnout_x = None
    burnout_z = None
    burnout_vx = None
    burnout_vz = None
    burnout_mass = None

    apogee_time = None
    apogee_x = None
    apogee_z = None

    chute_deploy_time = None
    chute_deploy_x = None
    chute_deploy_z = None
    parachute_deployed = False

    landing_time = None
    landing_x = None
    landing_speed = None
    landing_vx = None
    landing_vz = None

    max_alt = 0.0
    max_speed = 0.0
    max_mach = 0.0
    max_q = 0.0
    max_q_time = None
    max_q_alt = None
    max_cd = 0.0
    max_alpha_deg = 0.0

    prev_vz = vz

    ### Main time loop ###
    while t < t_max:
        ### Current state values ###
        current_x = state[0]
        current_z = max(0.0, state[1])
        current_vx = state[2]
        current_vz = state[3]
        current_theta = state[4]
        current_omega = state[5]
        current_mass = state[6]

        speed_sq = current_vx**2 + current_vz**2
        current_speed = sqrt(speed_sq)

        ### Aero values ###
        aero_x, aero_z, moment, current_mach, current_cd, current_q, \
            current_alpha, _ = aero_forces_and_moment(
                current_z,
                current_vx,
                current_vz,
                current_theta,
                vehicle.rocket,
                parachute,
                environment,
                atmosphere_func,
                parachute_deployed,
            )

        gravity_z = signed_gravity_force(current_mass, current_z)

        ### Thrust ###
        if parachute_deployed:
            thrust_x = 0.0
            thrust_z = 0.0
        else:
            thrust_mag = signed_thrust_force(t, vehicle.motor)
            thrust_x = thrust_mag * cos(current_theta)
            thrust_z = thrust_mag * sin(current_theta)

        ### Current acceleration ###
        force_x = thrust_x + aero_x
        force_z = thrust_z + aero_z + gravity_z

        ax = force_x / current_mass
        az = force_z / current_mass

        accel_sq = ax**2 + az**2
        current_accel = sqrt(accel_sq)

        ### Store histories ###
        time_hist.append(t)
        x_hist.append(current_x)
        z_hist.append(current_z)
        vx_hist.append(current_vx)
        vz_hist.append(current_vz)
        speed_hist.append(current_speed)
        theta_hist.append(degrees(current_theta))
        omega_hist.append(degrees(current_omega))
        alpha_hist.append(degrees(current_alpha))
        mass_hist.append(current_mass)
        accel_hist.append(current_accel)
        mach_hist.append(current_mach)
        q_hist.append(current_q)
        cd_hist.append(current_cd)

        ### Max metrics ###
        if current_z > max_alt:
            max_alt = current_z

        if current_speed > max_speed:
            max_speed = current_speed

        if current_mach > max_mach:
            max_mach = current_mach

        if current_q > max_q:
            max_q = current_q
            max_q_time = t
            max_q_alt = current_z

        if current_cd > max_cd:
            max_cd = current_cd

        alpha_deg = abs(degrees(current_alpha))
        if alpha_deg > max_alpha_deg:
            max_alpha_deg = alpha_deg

        ### Burnout ###
        if burnout_time is None and t >= vehicle.motor.burn_time_s:
            burnout_time = t
            burnout_x = current_x
            burnout_z = current_z
            burnout_vx = current_vx
            burnout_vz = current_vz
            burnout_mass = current_mass

        ### Apogee ###
        if prev_vz > 0.0 and current_vz <= 0.0 and apogee_time is None:
            apogee_time = t
            apogee_x = current_x
            apogee_z = current_z

        ### Parachute deployment ###
        if apogee_time is not None and not parachute_deployed \
            and t >= apogee_time + parachute.deploy_delay_s:
            parachute_deployed = True
            chute_deploy_time = t
            chute_deploy_x = current_x
            chute_deploy_z = current_z

        prev_vz = current_vz

        ### Save pre-step values ###
        prev_state = state.copy()
        prev_t = t

        ### RK4 integration ###
        state = rk4_step(
            state_derivatives,
            t,
            state,
            dt,
            vehicle,
            launch,
            parachute,
            environment,
            atmosphere_func,
            parachute_deployed,
        )

        ### Physical clamps ###
        state[1] = max(0.0, state[1])

        dry_mass = vehicle.rocket.dry_mass_kg
        if state[6] < dry_mass:
            state[6] = dry_mass

        t = t + dt

        ### Landing interpolation ###
        if t > 0.0 and prev_state[1] > 0.0 and state[1] <= 0.0:
            z0 = prev_state[1]
            z1 = state[1]

            if z0 != z1:
                frac = z0 / (z0 - z1)
            else:
                frac = 0.0

            landing_time = prev_t + frac * dt
            landing_x = prev_state[0] + frac * (state[0] - prev_state[0])
            landing_vx = prev_state[2] + frac * (state[2] - prev_state[2])
            landing_vz = prev_state[3] + frac * (state[3] - prev_state[3])

            landing_speed_sq = landing_vx**2 + landing_vz**2
            landing_speed = sqrt(landing_speed_sq)
            break

    ### Return dictionary ###
    results = {
        "time": time_hist,
        "x": x_hist,
        "altitude": z_hist,
        "vx": vx_hist,
        "vz": vz_hist,
        "speed": speed_hist,
        "theta_deg": theta_hist,
        "omega_deg_s": omega_hist,
        "alpha_deg": alpha_hist,
        "mass": mass_hist,
        "acceleration": accel_hist,
        "mach": mach_hist,
        "dynamic_pressure": q_hist,
        "cd": cd_hist,
        "burnout_time": burnout_time,
        "burnout_x": burnout_x,
        "burnout_altitude": burnout_z,
        "burnout_vx": burnout_vx,
        "burnout_vz": burnout_vz,
        "burnout_mass": burnout_mass,
        "apogee_time": apogee_time,
        "apogee_x": apogee_x,
        "apogee_altitude": apogee_z,
        "chute_deploy_time": chute_deploy_time,
        "chute_deploy_x": chute_deploy_x,
        "chute_deploy_altitude": chute_deploy_z,
        "landing_time": landing_time,
        "landing_x": landing_x,
        "landing_vx": landing_vx,
        "landing_vz": landing_vz,
        "landing_speed": landing_speed,
        "max_altitude": max_alt,
        "max_speed": max_speed,
        "max_mach": max_mach,
        "max_dynamic_pressure": max_q,
        "max_q_time": max_q_time,
        "max_q_altitude": max_q_alt,
        "max_cd": max_cd,
        "max_alpha_deg": max_alpha_deg,
    }

    return results

###########################################################################################
def run_monte_carlo(vehicle: Vehicle, launch: Launch, parachute: Parachute,
    environment: Environment, atmosphere_func, n_trials=200, dt=0.01,
    t_max=300.0):
    results_list = []

    for i in range(n_trials):
        ### Sample perturbations ###
        angle_perturb = random.uniform(-1.0, 1.0)
        launch_angle = launch.angle_deg + angle_perturb

        mass_scale = random.uniform(0.97, 1.03)
        dry_mass = vehicle.rocket.dry_mass_kg * mass_scale

        thrust_scale = random.uniform(0.95, 1.05)
        thrust = vehicle.motor.thrust_N * thrust_scale

        wind_scale = random.uniform(0.8, 1.2)

        delay_perturb = random.uniform(-0.5, 0.5)
        deploy_delay = parachute.deploy_delay_s + delay_perturb
        deploy_delay = max(0.0, deploy_delay)

        ### Build perturbed objects ###
        trial_rocket = Rocket(
            dry_mass_kg=dry_mass,
            diameter_m=vehicle.rocket.diameter_m,
            length_m=vehicle.rocket.length_m,
            mach_points=vehicle.rocket.mach_points[:],
            cd_points=vehicle.rocket.cd_points[:],
            cn_alpha=vehicle.rocket.cn_alpha,
            cm_alpha=vehicle.rocket.cm_alpha,
            moi_kg_m2=vehicle.rocket.moi_kg_m2,
        )

        trial_motor = Motor(
            thrust_N=thrust,
            burn_time_s=vehicle.motor.burn_time_s,
            prop_mass_kg=vehicle.motor.prop_mass_kg,
        )

        trial_vehicle = Vehicle(rocket=trial_rocket, motor=trial_motor)
        trial_launch = Launch(angle_deg=launch_angle)

        trial_parachute = Parachute(
            cd=parachute.cd,
            area_m2=parachute.area_m2,
            deploy_delay_s=deploy_delay,
            rot_damping_1_s=parachute.rot_damping_1_s,
            rot_stiffness_1_s2=parachute.rot_stiffness_1_s2,
            theta_eq_rad=parachute.theta_eq_rad,
        )

        trial_environment = scaled_environment(environment, wind_scale)

        ### Run one trial ###
        sim_result = simulate(
            trial_vehicle,
            trial_launch,
            trial_parachute,
            trial_environment,
            atmosphere_func,
            dt=dt,
            t_max=t_max,
        )

        ### Store selected results ###
        trial_result = {
            "max_altitude": sim_result["max_altitude"],
            "apogee_x": sim_result["apogee_x"],
            "landing_x": sim_result["landing_x"],
            "landing_speed": sim_result["landing_speed"],
            "max_dynamic_pressure": sim_result["max_dynamic_pressure"],
            "max_alpha_deg": sim_result["max_alpha_deg"],
            "launch_angle_deg": launch_angle,
            "dry_mass_kg": dry_mass,
            "thrust_N": thrust,
            "wind_scale": wind_scale,
            "deploy_delay_s": deploy_delay,
        }

        results_list.append(trial_result)

    return results_list

###########################################################################################
def print_monte_carlo_summary(mc_results):
    ### Gather valid results ###
    apogees = []
    landings = []
    landing_speeds = []

    for i in range(len(mc_results)):
        r = mc_results[i]

        if r["max_altitude"] is not None:
            apogees.append(r["max_altitude"])

        if r["landing_x"] is not None:
            landings.append(r["landing_x"])

        if r["landing_speed"] is not None:
            landing_speeds.append(r["landing_speed"])

    ### Output ###
    print("Monte Carlo Summary")
    print("-------------------")
    print(f"Trials: {len(mc_results)}")

    if apogees:
        print(f"Apogee mean: {statistics.mean(apogees):.2f} m")
        print(f"Apogee stdev: {statistics.pstdev(apogees):.2f} m")
        print(f"Apogee min: {min(apogees):.2f} m")
        print(f"Apogee max: {max(apogees):.2f} m")

    if landings:
        print(f"Landing x mean: {statistics.mean(landings):.2f} m")
        print(f"Landing x stdev: {statistics.pstdev(landings):.2f} m")
        print(f"Landing x min: {min(landings):.2f} m")
        print(f"Landing x max: {max(landings):.2f} m")

    if landing_speeds:
        print(f"Landing speed mean: {statistics.mean(landing_speeds):.2f} m/s")
        print(f"Landing speed stdev: {statistics.pstdev(landing_speeds):.2f} m/s")

###########################################################################################
def plot_monte_carlo_results(mc_results):
    ### Gather valid results ###
    apogees = []
    landings = []
    landing_speeds = []

    for i in range(len(mc_results)):
        r = mc_results[i]

        if r["max_altitude"] is not None:
            apogees.append(r["max_altitude"])

        if r["landing_x"] is not None:
            landings.append(r["landing_x"])

        if r["landing_speed"] is not None:
            landing_speeds.append(r["landing_speed"])

    ### Landing dispersion ###
    plt.figure()
    plt.scatter(landings, landing_speeds)
    plt.xlabel("Landing Downrange [m]")
    plt.ylabel("Landing Speed [m/s]")
    plt.title("Monte Carlo Landing Dispersion")
    plt.grid(True)

    ### Apogee histogram ###
    plt.figure()
    plt.hist(apogees, bins=20)
    plt.xlabel("Apogee [m]")
    plt.ylabel("Count")
    plt.title("Monte Carlo Apogee Distribution")
    plt.grid(True)

    ### Landing histogram ###
    plt.figure()
    plt.hist(landings, bins=20)
    plt.xlabel("Landing Downrange [m]")
    plt.ylabel("Count")
    plt.title("Monte Carlo Landing Downrange Distribution")
    plt.grid(True)

    plt.show()

###########################################################################################
def plot_wind_profile(environment: Environment):
    ### Smooth profile points ###
    z_min = environment.wind_altitudes_m[0]
    z_max = environment.wind_altitudes_m[-1]

    z_samples = []
    wind_samples = []

    n_points = 200
    for i in range(n_points):
        z = z_min + (z_max - z_min) * i / (n_points - 1)
        w = wind_velocity_x(z, environment)

        z_samples.append(z)
        wind_samples.append(w)

    ### Plot ###
    plt.figure()
    plt.plot(wind_samples, z_samples, label="Wind Profile")
    plt.scatter(environment.wind_speeds_m_s, environment.wind_altitudes_m,
        label="Profile Points")

    plt.xlabel("Wind Speed [m/s]")
    plt.ylabel("Altitude [m]")
    plt.title("Wind Speed vs Altitude")
    plt.grid(True)
    plt.legend()
    plt.show()

###########################################################################################
def plot_results(results):
    ### Event times ###
    time = results["time"]
    burnout_t = results["burnout_time"]
    max_q_t = results["max_q_time"]
    deploy_t = results["chute_deploy_time"]
    landing_t = results["landing_time"]

    ### Trajectory ###
    plt.figure()
    plt.plot(results["x"], results["altitude"], label="Trajectory")

    if results["burnout_x"] is not None and results["burnout_altitude"] is not None:
        plt.scatter(results["burnout_x"], results["burnout_altitude"], label="Burnout")

    if results["apogee_x"] is not None and results["apogee_altitude"] is not None:
        plt.scatter(results["apogee_x"], results["apogee_altitude"], label="Apogee")

    if results["chute_deploy_x"] is not None and results["chute_deploy_altitude"] is not None:
        plt.scatter(results["chute_deploy_x"], results["chute_deploy_altitude"],
            label="Chute Deploy")

    if results["landing_x"] is not None:
        plt.scatter(results["landing_x"], 0.0, label="Landing")

    plt.xlabel("Downrange Distance [m]")
    plt.ylabel("Altitude [m]")
    plt.title("Trajectory")
    plt.grid(True)
    plt.legend()

    ### Altitude ###
    plt.figure()
    plt.plot(time, results["altitude"], label="Altitude")
    if burnout_t is not None:
        plt.axvline(burnout_t, linestyle="--", label="Burnout")
    if max_q_t is not None:
        plt.axvline(max_q_t, linestyle=":", label="Max-Q")
    if deploy_t is not None:
        plt.axvline(deploy_t, linestyle="-.", label="Chute Deploy")
    if landing_t is not None:
        plt.axvline(landing_t, linestyle="-", label="Landing")
    plt.xlabel("Time [s]")
    plt.ylabel("Altitude [m]")
    plt.title("Altitude vs Time")
    plt.grid(True)
    plt.legend()

    ### Horizontal velocity ###
    plt.figure()
    plt.plot(time, results["vx"], label="Vx")
    if burnout_t is not None:
        plt.axvline(burnout_t, linestyle="--", label="Burnout")
    if deploy_t is not None:
        plt.axvline(deploy_t, linestyle="-.", label="Chute Deploy")
    if landing_t is not None:
        plt.axvline(landing_t, linestyle="-", label="Landing")
    plt.xlabel("Time [s]")
    plt.ylabel("Horizontal Velocity [m/s]")
    plt.title("Horizontal Velocity vs Time")
    plt.grid(True)
    plt.legend()

    ### Vertical velocity ###
    plt.figure()
    plt.plot(time, results["vz"], label="Vz")
    if burnout_t is not None:
        plt.axvline(burnout_t, linestyle="--", label="Burnout")
    if deploy_t is not None:
        plt.axvline(deploy_t, linestyle="-.", label="Chute Deploy")
    if landing_t is not None:
        plt.axvline(landing_t, linestyle="-", label="Landing")
    plt.xlabel("Time [s]")
    plt.ylabel("Vertical Velocity [m/s]")
    plt.title("Vertical Velocity vs Time")
    plt.grid(True)
    plt.legend()

    ### Speed ###
    plt.figure()
    plt.plot(time, results["speed"], label="Speed")
    if burnout_t is not None:
        plt.axvline(burnout_t, linestyle="--", label="Burnout")
    if max_q_t is not None:
        plt.axvline(max_q_t, linestyle=":", label="Max-Q")
    if deploy_t is not None:
        plt.axvline(deploy_t, linestyle="-.", label="Chute Deploy")
    if landing_t is not None:
        plt.axvline(landing_t, linestyle="-", label="Landing")
    plt.xlabel("Time [s]")
    plt.ylabel("Speed [m/s]")
    plt.title("Speed vs Time")
    plt.grid(True)
    plt.legend()

    ### Pitch angle ###
    plt.figure()
    plt.plot(time, results["theta_deg"], label="Theta")
    if burnout_t is not None:
        plt.axvline(burnout_t, linestyle="--", label="Burnout")
    if deploy_t is not None:
        plt.axvline(deploy_t, linestyle="-.", label="Chute Deploy")
    plt.xlabel("Time [s]")
    plt.ylabel("Pitch Angle [deg]")
    plt.title("Pitch Angle vs Time")
    plt.grid(True)
    plt.legend()

    ### Pitch rate ###
    plt.figure()
    plt.plot(time, results["omega_deg_s"], label="Omega")
    if burnout_t is not None:
        plt.axvline(burnout_t, linestyle="--", label="Burnout")
    if deploy_t is not None:
        plt.axvline(deploy_t, linestyle="-.", label="Chute Deploy")
    plt.xlabel("Time [s]")
    plt.ylabel("Pitch Rate [deg/s]")
    plt.title("Pitch Rate vs Time")
    plt.grid(True)
    plt.legend()

    ### Angle of attack ###
    plt.figure()
    plt.plot(time, results["alpha_deg"], label="Alpha")
    if burnout_t is not None:
        plt.axvline(burnout_t, linestyle="--", label="Burnout")
    if max_q_t is not None:
        plt.axvline(max_q_t, linestyle=":", label="Max-Q")
    if deploy_t is not None:
        plt.axvline(deploy_t, linestyle="-.", label="Chute Deploy")
    plt.xlabel("Time [s]")
    plt.ylabel("Angle of Attack [deg]")
    plt.title("Angle of Attack vs Time")
    plt.grid(True)
    plt.legend()

    ### Acceleration ###
    plt.figure()
    plt.plot(time, results["acceleration"], label="Acceleration")
    if burnout_t is not None:
        plt.axvline(burnout_t, linestyle="--", label="Burnout")
    if deploy_t is not None:
        plt.axvline(deploy_t, linestyle="-.", label="Chute Deploy")
    if landing_t is not None:
        plt.axvline(landing_t, linestyle="-", label="Landing")
    plt.xlabel("Time [s]")
    plt.ylabel("Acceleration [m/s²]")
    plt.title("Acceleration vs Time")
    plt.grid(True)
    plt.legend()

    ### Mass ###
    plt.figure()
    plt.plot(time, results["mass"], label="Mass")
    if burnout_t is not None:
        plt.axvline(burnout_t, linestyle="--", label="Burnout")
    plt.xlabel("Time [s]")
    plt.ylabel("Mass [kg]")
    plt.title("Mass vs Time")
    plt.grid(True)
    plt.legend()

    ### Mach ###
    plt.figure()
    plt.plot(time, results["mach"], label="Mach")
    if burnout_t is not None:
        plt.axvline(burnout_t, linestyle="--", label="Burnout")
    if max_q_t is not None:
        plt.axvline(max_q_t, linestyle=":", label="Max-Q")
    if deploy_t is not None:
        plt.axvline(deploy_t, linestyle="-.", label="Chute Deploy")
    plt.xlabel("Time [s]")
    plt.ylabel("Mach Number [-]")
    plt.title("Mach vs Time")
    plt.grid(True)
    plt.legend()

    ### Dynamic pressure ###
    plt.figure()
    plt.plot(time, results["dynamic_pressure"], label="Dynamic Pressure")
    if max_q_t is not None:
        plt.axvline(max_q_t, linestyle=":", label="Max-Q")
        plt.scatter(results["max_q_time"], results["max_dynamic_pressure"],
            label="Max-Q Point")
    if burnout_t is not None:
        plt.axvline(burnout_t, linestyle="--", label="Burnout")
    if deploy_t is not None:
        plt.axvline(deploy_t, linestyle="-.", label="Chute Deploy")
    plt.xlabel("Time [s]")
    plt.ylabel("Dynamic Pressure [Pa]")
    plt.title("Dynamic Pressure vs Time")
    plt.grid(True)
    plt.legend()

    ### Cd ###
    plt.figure()
    plt.plot(time, results["cd"], label="Cd")
    if burnout_t is not None:
        plt.axvline(burnout_t, linestyle="--", label="Burnout")
    if max_q_t is not None:
        plt.axvline(max_q_t, linestyle=":", label="Max-Q")
    if deploy_t is not None:
        plt.axvline(deploy_t, linestyle="-.", label="Chute Deploy")
    plt.xlabel("Time [s]")
    plt.ylabel("Drag Coefficient [-]")
    plt.title("Cd vs Time")
    plt.grid(True)
    plt.legend()

    plt.show()

###########################################################################################
def print_summary(results):
    ### Burnout speed ###
    burnout_speed = None
    if results["burnout_vx"] is not None and results["burnout_vz"] is not None:
        burnout_speed_sq = results["burnout_vx"]**2 + results["burnout_vz"]**2
        burnout_speed = sqrt(burnout_speed_sq)

    ### Output ###
    print("Rocket Flight Summary")
    print("---------------------")

    if results["burnout_time"] is not None:
        print(f"Burnout time: {results['burnout_time']:.2f} s")
        print(f"Burnout downrange: {results['burnout_x']:.2f} m")
        print(f"Burnout altitude: {results['burnout_altitude']:.2f} m")
        print(f"Burnout Vx: {results['burnout_vx']:.2f} m/s")
        print(f"Burnout Vz: {results['burnout_vz']:.2f} m/s")
        if burnout_speed is not None:
            print(f"Burnout speed: {burnout_speed:.2f} m/s")
        print(f"Burnout mass: {results['burnout_mass']:.2f} kg")

    if results["apogee_time"] is not None:
        print(f"Apogee time: {results['apogee_time']:.2f} s")
        print(f"Apogee downrange: {results['apogee_x']:.2f} m")
        print(f"Apogee altitude: {results['apogee_altitude']:.2f} m")

    if results["chute_deploy_time"] is not None:
        print(f"Chute deploy time: {results['chute_deploy_time']:.2f} s")
        print(f"Chute deploy downrange: {results['chute_deploy_x']:.2f} m")
        print(f"Chute deploy altitude: {results['chute_deploy_altitude']:.2f} m")

    if results["landing_time"] is not None:
        print(f"Landing time: {results['landing_time']:.2f} s")
        print(f"Landing downrange: {results['landing_x']:.2f} m")
        if results["landing_vx"] is not None:
            print(f"Landing Vx: {results['landing_vx']:.2f} m/s")
        if results["landing_vz"] is not None:
            print(f"Landing Vz: {results['landing_vz']:.2f} m/s")
        print(f"Landing speed: {results['landing_speed']:.2f} m/s")

    print(f"Max altitude: {results['max_altitude']:.2f} m")
    print(f"Max speed: {results['max_speed']:.2f} m/s")
    print(f"Max Mach: {results['max_mach']:.3f}")
    print(f"Max dynamic pressure: {results['max_dynamic_pressure']:.2f} Pa")

    if results["max_q_time"] is not None:
        print(f"Max-Q time: {results['max_q_time']:.2f} s")
        print(f"Max-Q altitude: {results['max_q_altitude']:.2f} m")

    print(f"Max effective Cd used in sim: {results['max_cd']:.3f}")
    print(f"Max |alpha|: {results['max_alpha_deg']:.2f} deg")

###########################################################################################
if __name__ == "__main__":
    ### Vehicle ###
    rocket = Rocket(
        dry_mass_kg=12.0,
        diameter_m=0.15,
        length_m=2.0,
        mach_points=[0.0, 0.3, 0.8, 1.0, 1.2, 2.0],
        cd_points=[0.45, 0.42, 0.48, 0.62, 0.55, 0.40],
        cn_alpha=2.0,
        cm_alpha=-0.8,
        moi_kg_m2=8.0
    )

    motor = Motor(thrust_N=850.0, burn_time_s=3.2, prop_mass_kg=4.0)
    vehicle = Vehicle(rocket=rocket, motor=motor)

    ### Launch / recovery / environment ###
    launch = Launch(angle_deg=85.0)

    parachute = Parachute(
        cd=1.5,
        area_m2=0.8,
        deploy_delay_s=1.0,
        rot_damping_1_s=2.0,
        rot_stiffness_1_s2=1.5,
        theta_eq_rad=-pi / 2
    )

    environment = Environment(
        wind_altitudes_m=[0.0, 100.0, 300.0, 600.0, 1000.0],
        wind_speeds_m_s=[2.0, 4.0, 6.0, 8.0, 10.0]
    )

    from atmosphere_model import atmosphere

    ### Nominal simulation ###
    results = simulate(vehicle, launch, parachute, environment, atmosphere,
        dt=0.01)

    ### Monte Carlo ###
    mc_results = run_monte_carlo(
        vehicle,
        launch,
        parachute,
        environment,
        atmosphere,
        n_trials=200,
        dt=0.01,
        t_max=300.0,
    )

    ### Output ###
    print_summary(results)
    plot_results(results)
    plot_wind_profile(environment)
    print_monte_carlo_summary(mc_results)
    plot_monte_carlo_results(mc_results)
