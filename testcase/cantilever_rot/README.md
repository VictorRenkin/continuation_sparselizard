## Data Source

The data used in this work originate from the following article:  
A. Martin et al., “Reduced-order Modeling of Geometrically Nonlinear Rotating Structures Using the Direct Parametrisation of Invariant Manifolds”, *Journal of Theoretical, Computational and Applied Mechanics*, 2023. DOI: 10.46298/jtcam.10430. -> Figure 3 (a) -> NNM using Plotdigitalizer.

## Simulation Setup

The harmonic order used in this simulation is set to 5.

### Physical Regions

| Physical Region                | Assigned IDs |
|--------------------------------|--------------|
| Physical region clamped        | 102          |
| Physical region material       | 104          |
| Physical region measured       | 103          |
| Physical region force          | 103          |

### Material Properties

| Property        | Value   | Unit  |
|-----------------|---------|-------|
| Young's modulus | 1.04e11 | Pa    |
| Poisson's ratio | 0.3     | –     |
| Density         | 4400    | kg/m³ |

Due to the centrifugal stiffening induced by rotation, the dynamic characteristics of the system evolve with the rotational speed. In order to maintain a comparable level of vibration amplitude in the nonlinear frequency response analyses, the forcing amplitude is increased as the rotational speed increases. In parallel, the damping properties are adjusted to ensure a consistent level of energy dissipation. Specifically, the damping ratio is increased such that it remains equal to $0.05\%$ for all considered rotational speeds.

### Damping and Forcing Parameters

| RPM  | Alpha damping                  | Applied force              |
|------|--------------------------------|----------------------------|
| 0    | 4 · π · 0.0005 · 15.7472        | (-2, 0, 0)                 |
| 500  | 4 · π · 0.0005 · 16.4969        | (-2, 0, 0) · 1.1027        |
| 1000 | 4 · π · 0.0005 · 18.4749        | (-2, 0, 0) · 1.397         |
| 1500 | 4 · π · 0.0005 · 21.1705        | (-2, 0, 0) · 8.89 / 4.77   |
| 2000 | 4 · π · 0.0005 · 24.2064        | (-2, 0, 0) · 11.79 / 4.77  |