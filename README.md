# Two-Body Problem Simulation

## Overview
This repository contains a MATLAB-based simulation of the *Two-Body Problem*, which models the motion of a satellite in orbit around Earth. The project implements both analytical and numerical methods to propagate satellite trajectories using Kepler's laws and Newtonian mechanics. Additionally, it explores attitude representations and torque-free motion using Euler's equations.

This project was developed as part of the coursework for the module *Space Dynamics and Missions* of the *MSc in Space Engineering* at the *University of Surrey*.

## Objectives
- Solve Kepler’s equation using Newton’s method to determine the Eccentric Anomaly.
- Propagate satellite motion using both analytical methods and numerical integration in the Earth-Centered Inertial (ECI) frame.
- Convert orbital elements to position and velocity vectors in the ECI frame.
- Integrate the equations of motion in the Earth-Centered Earth-Fixed (ECEF) frame.
- Model satellite attitude using Direction Cosine Matrices (DCM), Euler angles, and quaternions.
- Simulate torque-free motion using Euler’s equations.

## File Structure
```
📁 two-body-problem
│── 📂 src                       # MATLAB source files
│   │── Kepler.m                  # Solves Kepler's equation using Newton's method
│   │── COE2RV.m                  # Converts classical orbital elements to ECI position/velocity
│   │── TBP_ECI.m                 # Integrates the equations of motion in ECI frame
│   │── TBP_ECEF.m                # Integrates the equations of motion in ECEF frame
│   │── AttitudeDynamics.m        # Solves Euler's equations for torque-free motion
│   │── main.m                    # Main script for running the simulation
│   │── BlueMarble_square.png     # Image used for smiulations
│── README.md                     # Project documentation
```

## Getting Started
### Prerequisites
- MATLAB
- Basic knowledge of orbital mechanics and attitude dynamics

### Running the Simulation
1. Clone the repository:
   ```sh
   git clone https://github.com/JorgeACM/two-body-problem.git
   cd two-body-problem/src
   ```
2. Open MATLAB and navigate to the `src` directory.
3. Run the main script:
   ```matlab
   main
   ```
4. The script will generate orbital trajectory plots, attitude evolution, and numerical validation results.

## Features
### 1. **Solving Kepler’s Equation**
   - Uses Newton’s method to compute the Eccentric Anomaly.
   - Converts Mean Anomaly to True Anomaly.

### 2. **Orbit Propagation**
   - Analytical propagation using Keplerian elements.
   - Numerical integration in the ECI frame with `ode113`.
   - Comparison of analytical vs. numerical results.

### 3. **ECEF Frame Integration**
   - Converts initial conditions from ECI to ECEF.
   - Accounts for Earth’s rotation and Coriolis effects.

### 4. **Attitude Representations**
   - Computes Direction Cosine Matrices (DCM) for different reference frames.
   - Converts between Euler angles, quaternions, and principal axes.
   - Simulates attitude kinematics and torque-free motion.

### 5. **Torque-Free Motion Simulation**
   - Solves Euler’s equations for rotational dynamics.
   - Validates conservation of angular momentum and kinetic energy.

## Results
- 3D visualizations of satellite trajectories in both ECI and ECEF frames.
- Plots of Mean, Eccentric, and True Anomalies over time.
- Euler angles and angular velocity time histories.
- Energy conservation validation for torque-free motion.

## References
- Curtis, H. D. *Orbital Mechanics for Engineering Students*, 4th Edition.
- Vallado, D. A. *Fundamentals of Astrodynamics and Applications*.
- Schaub, H., Junkins, J. *Analytical Mechanics of Space Systems*.

## License
This project is licensed under the MIT License.

For more information, please refer to the report available in the GitHub repository.
