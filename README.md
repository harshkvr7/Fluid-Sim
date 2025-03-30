# 2D Fluid Simulation with SDL2

This repository contains a simple 2D fluid simulation implemented in C with real-time visualization using SDL2. The simulation models fluid behavior by simulating the diffusion and advection of a "dye" (density) within a fluid field. Users can interact with the simulation by injecting density and velocity forces using the mouse.

## Overview

The project combines a fluid simulation algorithm with SDL2-based visualization. Key features include:

- **Fluid Simulation:**  
  The simulation is based on a grid where each cell stores density and velocity components. The algorithm updates the state using:
  - **Diffusion:** Models the spread of dye and velocity across the grid.
  - **Advection:** Moves dye and velocity along the fluid flow.
  - **Projection:** Enforces the incompressibility of the fluid by ensuring the velocity field is divergence-free.

- **SDL2 Visualization:**  
  The simulation is rendered in a window where each grid cell is represented by a colored rectangle. The color gradient (from blue for low density to red for high density) visually represents the concentration of dye. Users can interact with the simulation by clicking and dragging the mouse to add density and force.

## Features

- Real-time 2D fluid simulation.
- Interactive dye and force injection using the mouse.
- Visual representation of fluid density using a blue-to-red color gradient.
- Parameter tuning for grid size, diffusion, viscosity, and time step.

## Requirements

- [SDL2](https://www.libsdl.org/) library installed on your system.
- A C compiler (e.g., GCC) with math library support.

## Installation

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/harshkvr7/Fluid-Sim.git
   cd Fluid-sim

2. **Compile and Run:**

   ```bash
   gcc -o fluid fluid.c -lSDL2 -lm && ./fluid

## Reference 

This project is inspired by the article [Fluid Simulation for Dummies](https://www.mikeash.com/pyblog/fluid-simulation-for-dummies.html) by Mike Ash. The article provides an accessible introduction to fluid simulation techniques and was a major reference in implementing the simulation algorithm in this project.
