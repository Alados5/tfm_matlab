# Robotic Cloth Manipulation: Real Implementation with MPC and RL
**Double Master’s Thesis - ETSEIB, 2021**

This repository includes the complete source code written in Matlab used to execute closed-loop simulations of a robot moving a cloth piece using MPC, and also to improve this controller with Reinforcement Learning techniques.

**Affiliation**: [Institut de Robòtica i Informàtica Industrial (IRI) CSIC-UPC](https://www.iri.upc.edu/), Barcelona.

## Main Structure
- **closed_loop_sims**: Main codes to execute closed-loop simulations
- **data**: Some saved results and plotting codes, and also the files of the closed-loop reference trajectories
- **learn_control**: Codes and results obtained applying Reinforcement Learning to tune the MPC
- **learn_model**: Codes and results obtained applying Reinforcement Learning to find the parameters of the cloth model
- **mpc_progression**: Snippets of code representing the progressive changes done during the project, mentioned in its Memory
- **required_files**: CasADi toolbox, needed for the optimizations, and additional functions related to the cloth models
- **with_robot**: Codes using the Robotics Toolbox to initialize and plot a WAM robot in Matlab
