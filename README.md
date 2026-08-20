# PoincareMaps.jl

Numerical computation of **Poincaré maps** for dynamical systems using TaylorIntegration.jl

---

## Overview

`PoincareMaps.jl` is a Julia package for computing intersections of trajectories with a Poincaré surface of section.

The package is designed around Taylor-series integration and is particularly aimed at Hamiltonian dynamical systems, where Poincaré sections are a useful tool for studying periodic orbits, invariant structures, resonances, and chaotic dynamics.

The package currently builds on [`TaylorIntegration.jl`](https://github.com/Perez-Hernandez/TaylorIntegration.jl) and [`TaylorSeries.jl`](https://github.com/JuliaDiff/TaylorSeries.jl).

---

## Main Features

- **Taylor-series integration** of dynamical systems.
- **Detection of crossings** with a user-defined surface.
- **Optional additional conditions** for accepting a crossing.
- **Newton refinement** of the crossing time.
- **User-defined boundary-condition transformations** after each integration step.
- **Optional stopping conditions** during integration.
- **Energy-constrained generation** of initial conditions.
- **Computation of Poincaré sections** from multiple initial conditions.
- Support for **arbitrary state dimension**.
- Designed to work with Hamiltonian systems, but **not restricted to them**.

---

## Mathematical Formulation

Consider an autonomous dynamical system:

\[ \dot{x} = f(x), \quad x \in \mathbb{R}^n \]

A Poincaré section is defined through a scalar function:

$$g(x) = 0$$

During integration, the package detects changes in the sign of $g$ between consecutive Taylor integration steps. Once a crossing is detected, its location is refined using Newton-Raphson iterations applied to the Taylor polynomial.

### Surface Function Interface

The surface function has the form:

```julia
function g(dx, x, params, t)
    return (condition, g_value)
end
```

Where:
- `condition::Bool`: Determines whether the crossing satisfies additional user-defined requirements.
- `g_value`: The quantity whose zero defines the surface.

#### Examples

1. **Simple Section at $q = 0$:**
   ```julia
   function g(dx, x, params, t)
       return (true, x[1])
   end
   ```

2. **Direction-Restricted Section ($q = 0$ with $\dot{q} > 0$):**
   ```julia
   function g(dx, x, params, t)
       condition = dx[1] > 0
       return (condition, x[1])
   end
   ```

> The second component is represented internally by a Taylor polynomial, allowing the crossing time to be refined without performing a separate integration.

---

## Basic Usage

A minimal example for the harmonic oscillator ($\dot{q} = p$, $\dot{p} = -q$):

```julia
using PoincareMaps

function f!(dx, x, params, t)
    dx[1] = x[2]
    dx[2] = -x[1]
end

function bc!(x, params, t)
    nothing
end

function g(dx, x, params, t)
    return (true, x[1])
end

q0 = [1.0, 0.0]
t0 = 0.0
tmax = 10π

abstol = 1e-16
order = 20
maxevents = 20
maxsteps = 1000
params = nothing

cache = PoincareMaps.init_cache_ps(
    t0,
    q0,
    maxevents,
    order,
    f!,
    params;
    parse_eqs=true,
)

nevents = PoincareMaps.taylorinteg_ps!(
    f!,
    bc!,
    g,
    q0,
    t0,
    tmax,
    abstol,
    cache,
    params;
    maxsteps=maxsteps,
    maxevents=maxevents,
    eventorder=0,
)

crossings = cache.xv[:, 1:nevents-1]
```

For this system and section, the intersections alternate approximately between:
- $(q, p) = (0, -1)$
- $(q, p) = (0, +1)$

*Note: The implementation uses an internal event counter that starts at 1, so the number of detected crossings is `nevents - 1`.*

---

## PoincareMap

The higher-level `PoincareMap` function allows a collection of initial conditions to be generated and evaluated on a fixed energy surface.

### Interface

```julia
PoincareMap(
    f!,
    bc!,
    g,
    H,
    nrfind,
    seed,
    E,
    params,
    tmax,
    abstol,
    order;
    parse_eqs=true,
    maxsteps=500,
    maxevents=500,
    eventorder=0,
    newtoniter=10,
    nrabstol=eps(E),
    Etol=eps(E),
)
```

### Parameter Description

| Parameter | Description |
| :--- | :--- |
| `f!` | The vector field function |
| `bc!` | Applies boundary-condition or coordinate-normalization transformations |
| `g` | Defines the Poincaré section and crossing conditions |
| `H` | The Hamiltonian or conserved quantity used to impose target energy |
| `nrfind` | Determines coordinate & interval used to construct initial conditions on energy surface |
| `seed` | Contains the initial-condition seeds |
| `E` | Target energy |
| `params` | Model parameters |
| `tmax` | Maximum integration time |
| `abstol` | Taylor integration absolute tolerance |
| `order` | Taylor expansion order |

The returned object is a matrix whose rows correspond to detected intersections:

$$	ext{number\_of\_crossings} 	imes 	ext{state\_dimension}$$

Thus, for a four-dimensional Hamiltonian system, every row represents one point of the Poincaré section in the full phase space.

---

## Energy-Constrained Initial Conditions

For Hamiltonian systems, `PoincareMap` can generate initial conditions satisfying $H(x) = E$.

The function `nrfind` provides the coordinate and admissible interval over which the initial condition is adjusted.

```julia
function nrfind(x, params, t)
    xmin = 0.0
    xmax = 1.0
    Δx = (xmax - xmin) / 50

    return (4, Δx, xmax, xmin)
end
```

The selected coordinate is then refined with Newton's method until the Hamiltonian satisfies the requested energy within `Etol`. This allows a set of seeds to be projected onto a specified energy surface before trajectory integration.

---

## Poincaré Section Conditions

The Boolean returned by `g` can be used to distinguish geometrically different crossings.

### Examples

1. **Directional Selection:**
   ```julia
   function g(dx, x, params, t)
       condition = dx[1] > 0
       return (condition, x[1])
   end
   ```
   *Detects the surface $x_1 = 0$, but retains only crossings in the positive $x_1$ direction.*

2. **Complex State-Dependent Selection:**
   ```julia
   function g(dx, x, params, t)
       condition = x[2] > 0 && abs(x[3]) < 1
       return (condition, x[1])
   end
   ```

This mechanism is useful for selecting a particular branch of a section or avoiding duplicate intersections.

---

## Boundary Conditions

The `bc!` function is applied to the state after every Taylor integration step:

```julia
function bc!(x, params, t)
    # modify x in place
    nothing
end
```

This is useful for systems with periodic coordinates (e.g., wrapping angles):

```julia
function bc!(x, params, t)
    x[1] = mod(x[1] + π, 2π) - π
    nothing
end
```

---

## Integration Architecture

The low-level integration functionality is separated into several components:

### `VectorCachePS`

Stores the Taylor integration state and the arrays required for event detection:

```julia
struct VectorCachePS{XV,XAUX,T,X,DX,RV,PARSE_EQS} <: AbstractVectorCache
    xv::XV
    xaux::XAUX
    t::T
    x::X
    dx::DX
    rv::RV
    parse_eqs::PARSE_EQS
end
```

Initialization:

```julia
cache = init_cache_ps(
    t0,
    q0,
    maxevents,
    order,
    f!,
    params;
    parse_eqs=true,
)
```

Detected crossings are stored in `cache.xv`. 

> **Note:** The package intentionally does not store the intersection times. The primary output of the Poincaré computation is the set of intersection points in phase space.

### `taylorinteg_ps!`

Performs integration and event detection with the following workflow:
1. Initialize the Taylor integration cache.
2. Generate the Taylor expansion of the vector field.
3. Advance the trajectory using `TaylorIntegration`.
4. Evaluate the surface function.
5. Detect crossings.
6. Refine crossing time using Newton-Raphson.
7. Apply boundary-condition function.
8. Store accepted intersections.

*Two variants are provided:*
- A standard integration variant.
- A variant accepting `lims`, which can terminate integration when a user-defined condition is reached.

### `findroot_ps!`

Performs crossing detection and Newton refinement. When a crossing is detected, linear interpolation provides an initial estimate for the crossing time. Newton-Raphson iterations then refine this estimate using the Taylor polynomial of the section function.

---

## Dependencies

- `TaylorIntegration.jl`
- `TaylorSeries.jl`
- `LinearAlgebra`

### Development Target
- `TaylorIntegration 0.18.3`
- Compatible `TaylorSeries` version (specified by package compatibility constraints).

---

## Installation

Once registered or available from its repository, install via Julia's package manager:

```julia
using Pkg
Pkg.add(url="https://github.com/LeosEquation/PoincareMaps.jl")
```

For local development:

```julia
using Pkg
Pkg.develop(path="path/to/PoincareMaps")
```

Load the package:

```julia
using PoincareMaps
```

---

## Tests

Uses Julia's standard `Test` framework with structure:

```text
test/
├── runtests.jl
├── integration/
│   └── taylorinteg_ps.jl
└── poincaremap/
    └── ...
```

Run test suite:

```julia
] test
```
*or*
```julia
Pkg.test()
```

The integration tests include a harmonic oscillator example to verify event detection and crossing accuracy against analytical solutions.

---

## Examples

Examples and exploratory notebooks are kept separately from the package source, demonstrating:
- Basic Poincaré sections
- Event detection
- Hamiltonian systems
- Energy-constrained initial conditions
- Periodic coordinates and boundary conditions
- Complex nonlinear dynamics

---

## Project Structure

```text
PoincareMaps.jl/
├── src/
│   ├── PoincareMaps.jl
│   ├── poincaremap.jl
│   └── integration/
│       ├── cache.jl
│       ├── findroot.jl
│       └── taylorinteg.jl
├── test/
│   ├── runtests.jl
│   ├── integration/
│   │   └── taylorinteg_ps.jl
│   └── poincaremap/
├── examples/
└── Project.toml
```

---

## Design Goals

- **Independence:** Keep Poincaré-map implementation independent from specific Hamiltonian models.
- **High Accuracy:** Reuse high-order Taylor integration machinery (`TaylorIntegration.jl`).
- **Refinement:** Refine intersections using Taylor polynomial information rather than low-order interpolation.
- **Composability:** Maintain a small and composable core interface.
- **Flexibility:** Allow user-defined energy constraints and section conditions.
- **Foundation:** Provide a framework for numerical studies of periodic orbits and bifurcations.

---

## Status

`PoincareMaps.jl` is under active development. Low-level integration is functional and tested; high-level energy-constrained Poincaré mapping features are being refined. The API may change prior to the first stable release.

---

## License & Author

- **Author:** Leonel Mayorga López ([GitHub](https://github.com/LeosEquation))
- **License:** Distributed under the license specified in the repository.
