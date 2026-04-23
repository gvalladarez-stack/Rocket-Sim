# Rocket Flight Dynamics Simulator
This project implements a 3DOF planar rocket flight dynamics simulator in Python. The model captures both translational motion and pitch dynamics using a 7-state system and propagates the trajectory using a Runge–Kutta 4 (RK4) integrator.

The simulation includes aerodynamic forces and moments based on Mach number and angle of attack, altitude-dependent atmospheric properties, and wind-relative flight dynamics. A recovery phase is modeled using a parachute with drag-based descent and damped rotational behavior.

In addition to single-flight simulation, the project includes a Monte Carlo analysis to evaluate sensitivity to uncertainty in key parameters such as launch angle, mass, thrust, wind, and parachute deployment timing. This allows estimation of landing dispersion, apogee variation, and overall flight variability.

This model is intended as a mid-fidelity engineering simulation that balances physical realism with computational simplicity. It is designed to demonstrate core concepts in flight dynamics, numerical integration, and uncertainty analysis, rather than to serve as a high-fidelity flight prediction tool.

Key Features:
- 3DOF planar rigid-body flight dynamics (translation + pitch rotation)
- RK4 numerical integration of a 7-state system
- Mach-dependent drag model with interpolated Cd data
- Angle-of-attack-based normal force and pitching moment modeling
- Altitude-dependent atmosphere and wind profile
- Parachute descent with drag and rotational damping/restoring behavior
- Event detection (burnout, max dynamic pressure, apogee, landing)
- Monte Carlo dispersion analysis for uncertainty quantification
# How to Run
Clone the repository and run the main script:

`python rocketsim.py`

The script will:
- run a nominal trajectory simulation
- generate trajectory and performance plots
- perform a Monte Carlo analysis and display dispersion results
# Notes
This project was developed to explore and demonstrate core aerospace modeling concepts, including the relationship between aerodynamic forces, vehicle stability, and environmental effects on flight trajectory. The emphasis is on clarity of implementation and physical reasoning rather than maximum model fidelity.
