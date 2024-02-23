``@meta
CurrentModule = GCIdentifier
```
# GCIdentifier.jl
An extensible [Julia](http://julialang.org) package for the modelling of fluids using thermodynamic equations of state. These include the standard cubics (van der Waals, Redlich-Kwong, Peng-Robinson, _etc._), SAFT-type equations (PC-SAFT, SAFT-VR Mie, SAFT-$\gamma$ Mie, _etc._), empirical equations (GERG2008, IAWPS95), Activity coefficient models (NRTL, UNIFAC, COSMO-SAC, _etc._) and many more.

The documentation is laid out as follows:

- **Group Assignment**: Find out how to assign groups to a species within a group-contribution method.
- **Finding Missing Groups**: Find out how to identify missing groups for a given species.
- **Custom Groups**: Find out how to implement your own groups within GCIdentifier.
- **API**: A list of all available methods.

### Authors

- [Pierre J. Walker](mailto:pjwalker@caltech.edu), California Institute of Technology
- [Andrés Riedemann](mailto:andres.riedemann@gmail.com), University of Concepción

### License

GCIdentifier.jl is licensed under the [MIT license](https://github.com/ClapeyronThermo/GCIdentifier.jl/blob/master/LICENSE.md).

### Installation

GCIdentifier.jl is a registered package, it can be installed from the general registry by:

```
pkg> add GCIdentifier
```