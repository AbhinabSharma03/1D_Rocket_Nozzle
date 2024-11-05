---

# 🚀 **FlowDynamics**: Quasi 1D Fluid Flow Through a Converging-Diverging Rocket Nozzle

**Description**  
This project performs Computational Fluid Dynamics (CFD) analysis of one-dimensional fluid flow through a converging-diverging rocket nozzle. It leverages MacCormack’s scheme to simulate both subsonic and supersonic conditions, using finite difference methods to solve the Euler equations for compressible flow.

**Key Features**:
- **Two Flow Regimes**: Simulates both subsonic and supersonic flow conditions.
- **Numerical Scheme**: Implements MacCormack's predictor-corrector scheme for improved stability and accuracy.
- **Results Visualization**: Plots pressure, density, temperature, and Mach number distributions along the nozzle length.

## 📁 Project Structure

- `Project_3.py` - Python script for the CFD simulation. Adjusts flow regimes based on the set conditions and provides plots for analysis.
- `project_3_cfd.pdf` - Detailed project report explaining methodology, equations, assumptions, and results.

## 🚀 Getting Started

### Prerequisites
- **Python 3.8+**
- **Libraries**: `numpy`, `matplotlib`

Install dependencies with:
```bash
pip install numpy matplotlib
```

### Usage
1. **Select Flow Regime**: Set `isentropic_subsonic = True` in `Project_3.py` for subsonic flow, or `False` for supersonic flow.
2. **Run Simulation**:
   ```bash
   python Project_3.py
   ```
3. **View Results**: The script generates plots for pressure, density, temperature, and Mach number distributions along the nozzle.

### Sample Output

- **Pressure Distribution** along the nozzle
- **Density and Temperature** profiles
- **Mach Number** changes, visualizing transition through sonic speeds.

## 📖 Methodology

The project uses MacCormack's scheme with the following features:
- **Predictor and Corrector Steps** to ensure accurate results in both flow regimes.
- **Area Distribution** based on the parabolic geometry function for converging-diverging nozzles.
- **Boundary and Initial Conditions** tailored for stable simulation across different cases.

## 📊 Results

The project provides comprehensive plots for each flow variable, allowing for detailed analysis of flow characteristics in different regimes. Validation against analytical results shows consistency, especially in capturing supersonic acceleration and subsonic recovery patterns.

## 📜 References
- Anderson, J. D., *Modern Compressible Flow: With Historical Perspective*

## 📈 Future Enhancements
- Integration of more complex nozzle geometries.
- Extension to three-dimensional CFD analysis.

--- 
