# Rocket Flight Dynamics Simulator
This project implements a 3DOF planar rocket flight dynamics simulator in Python. The model captures both translational motion and pitch dynamics using a 7-state system and propagates the trajectory using a Runge–Kutta 4 (RK4) integrator.

The simulation includes aerodynamic forces and moments based on Mach number and angle of attack, altitude-dependent atmospheric properties, and wind-relative flight dynamics. A recovery phase is modeled using a parachute with drag-based descent and damped rotational behavior.

In addition to single-flight simulation, the project includes a Monte Carlo analysis to evaluate sensitivity to uncertainty in key parameters such as launch angle, mass, thrust, wind, and parachute deployment timing. This allows estimation of landing dispersion, apogee variation, and overall flight variability.

This model is intended as a mid-fidelity engineering simulation that balances physical realism with computational simplicity. It is designed to demonstrate core concepts in flight dynamics, numerical integration, and uncertainty analysis, rather than to serve as a high-fidelity flight prediction tool.

## Key Features:
- 3DOF planar rigid-body flight dynamics (translation + pitch rotation)
- RK4 numerical integration of a 7-state system
- Mach-dependent drag model with interpolated C<sub>d</sub> data
- Angle-of-attack-based normal force and pitching moment modeling
- Altitude-dependent atmosphere and wind profile
- Parachute descent with drag and rotational damping/restoring behavior
- Event detection (burnout, max dynamic pressure, apogee, landing)
- Monte Carlo dispersion analysis for uncertainty quantification
## How to Run
Clone the repository and run the main script:

`python rocketsim.py`

The script will:
- run a nominal trajectory simulation
- generate trajectory and performance plots
- perform a Monte Carlo analysis and display dispersion results
## Notes
This project was developed to explore and demonstrate core aerospace modeling concepts, including the relationship between aerodynamic forces, vehicle stability, and environmental effects on flight trajectory. The emphasis is on clarity of implementation and physical reasoning rather than maximum model fidelity.
## Model Description & Limitations
Overview
This project implements a 3DOF planar rocket flight dynamics simulator with atmospheric modeling, aerodynamic forces and moments, and recovery system simulation. The model propagates a 7-state system [x, z, vx, vz, θ, ω, m], where translational motion and pitch dynamics are coupled through aerodynamic forces and moments. Numerical integration is performed using a fixed-step Runge–Kutta 4 (RK4) method.

1. Model Assumptions

Geometry & Motion:
Motion is restricted to a 2D vertical plane (x–z). Only pitch rotation is modeled (no yaw or roll). The rocket is treated as a rigid body with constant moment of inertia. The Earth is modeled as flat and non-rotating.

Atmosphere:
The simulation uses a standard atmosphere model (ISA-like). It assumes dry air, no humidity effects, and no temperature anomalies. Atmospheric properties vary only with altitude.

Wind:
Wind is modeled as horizontal only (x-direction), altitude-dependent, and time-invariant. The wind profile is defined by a user-provided table and linearly interpolated.

Propulsion:
Thrust is constant during burn and zero after burnout. The thrust vector is perfectly aligned with the rocket body axis. No thrust misalignment or gimbal effects are modeled.

Mass & Inertia:
Mass decreases linearly during burn and remains constant after burnout. The moment of inertia is fixed. No propellant slosh or center-of-mass shift is modeled.

2. Aerodynamics

Drag:
Drag is modeled using the standard equation D = 0.5 * ρ * V² * C<sub>d</sub> * A. The drag coefficient C<sub>d</sub> is Mach-dependent and defined via an interpolation table. Reynolds number effects and surface roughness are not modeled.

Normal Force (Angle of Attack):
The normal force coefficient is modeled as C<sub>n</sub> = C<sub>n<sub>α</sub></sub> · α. This is a linear aerodynamic model valid for small angles of attack.

Pitching Moment:
The pitching moment coefficient is modeled as C<sub>m</sub> = C<sub>m<sub>α</sub></sub> · α. This provides restoring torque for stable configurations. Coefficients are constant and do not vary with Mach or Reynolds number.

Angle of Attack:
Angle of attack is defined as alpha = theta − gamma, where gamma is the flight-path angle derived from air-relative velocity. This assumes attached flow and no flow separation.

Parachute Model:
Parachute deployment is instantaneous and triggered after apogee plus a specified delay. The parachute drag coefficient is constant during descent. No canopy inflation dynamics or oscillations are modeled.

Rotational Behavior Under Parachute:
Post-deployment rotational dynamics are modeled as a damped restoring system: θ¨ = −k(θ − θ_eq) − cω. This approximates stabilization under canopy but does not explicitly model pendulum motion or line dynamics.

Numerical Integration:
A fixed timestep RK4 integrator is used. Event detection (burnout, apogee, landing) is timestep-based. Landing time is refined using linear interpolation between steps.

3. Model Limitations

Flight Dynamics:
The model does not include yaw or roll dynamics and is not a full 6DOF simulation. Aerodynamic damping terms such as Cm_q are not included. There is no coupling between angular rate and aerodynamic force coefficients.

Aerodynamics:
The linear aerodynamic model becomes inaccurate at large angles of attack or in strongly transonic or shock-dominated regimes. Flow separation and nonlinear aerodynamic effects are not modeled. There is no explicit modeling of lift surfaces or fin geometry.

Propulsion:
The model assumes constant thrust and does not include a thrust curve. Ignition transients and thrust oscillations are not modeled.

Wind & Environment:
The model does not include turbulence, gusts, or time-varying wind. Only horizontal wind is considered, and no vertical wind component is included.

Recovery System:
Parachute inflation lag is not modeled. The model does not include dynamic canopy behavior, swinging or pendulum motion, line twist, or instability effects.

Ground Interaction:
The model does not include ground impact physics or terrain modeling.

 4. Monte-Carlo Analysis

Parameters Randomized:
Launch angle, rocket mass, motor thrust, wind profile scaling, and parachute deployment delay are varied within predefined ranges.

Outputs Analyzed:
Apogee altitude, landing downrange distance, landing speed, maximum dynamic pressure (Max-Q), and maximum angle of attack are recorded and analyzed.

Purpose:
The Monte Carlo analysis provides insight into flight variability, sensitivity to modeling assumptions, and a more realistic estimate of expected flight behavior.

