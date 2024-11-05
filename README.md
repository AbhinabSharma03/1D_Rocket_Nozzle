Here is the README in the requested GitHub format:

# 🚀 FlowDynamics: Quasi 1D Fluid Flow Through a Converging-Diverging Rocket Nozzle

This project performs a Computational Fluid Dynamics (CFD) analysis of one-dimensional fluid flow through a converging-diverging rocket nozzle. Using MacCormack's scheme and the finite difference method, it simulates both subsonic and supersonic flow conditions to solve the Euler equations for compressible flow.

## Features

- **Flow Regimes**: Simulates both isentropic subsonic and subsonic-supersonic flow conditions.
- **Numerical Method**: Implements MacCormack's predictor-corrector scheme for stable, accurate results.
- **Visualization**: Plots pressure, density, temperature, and Mach number distributions along the nozzle length.

## Project Structure

- **`Project_3.py`** - Python script for the CFD simulation with adjustable flow regimes and plots for analysis.
- **`project_3_cfd.pdf`** - Project report with a detailed explanation of methodology, equations, assumptions, and results.

## Getting Started

### Prerequisites

- **Python 3.8+**
- **Libraries**: `numpy`, `matplotlib`

Install dependencies by running:

```bash
pip install numpy matplotlib
```

### Usage

1. **Set Flow Regime**: Open `Project_3.py` and set `isentropic_subsonic = True` for subsonic flow or `False` for supersonic flow.
2. **Run Simulation**:

   ```bash
   python Project_3.py
   ```

3. **View Results**: The script generates plots for pressure, density, temperature, and Mach number along the nozzle.

### Example Output

- **Pressure Distribution** along the nozzle length.
- **Density and Temperature** profiles for visual analysis.
- **Mach Number** profiles indicating flow behavior through the nozzle.

## Methodology

This project uses MacCormack's scheme with:

- **Predictor and Corrector Steps** for accuracy in different flow regimes.
- **Variable Area Distribution** based on the parabolic geometry function for a converging-diverging nozzle.
- **Boundary Conditions** designed for simulation stability across cases.

## Results

The project generates comprehensive plots for each variable, allowing detailed flow analysis. The simulation results align well with analytical expectations, especially in supersonic acceleration and subsonic recovery patterns.

## References

- Anderson, J. D., *Modern Compressible Flow: With Historical Perspective*

## Future Enhancements

- Adding support for more complex nozzle geometries.
- Extending the simulation to three-dimensional CFD analysis.
