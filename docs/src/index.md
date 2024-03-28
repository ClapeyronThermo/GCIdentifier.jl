```@meta
CurrentModule = GCIdentifier
```
# GCIdentifier.jl
Welcome to GCIdentifier! This module provides utilities needed to fragment a given molecular SMILES (or name) based on the groups provided in existing group-contribution methods (such as UNIFAC, Joback's method and SAFT-$\gamma$ Mie). Additional functionalities have been provided to automatically identify and propose new groups.

Group-contribution approaches are vital when it comes to computer-aided molecular design (CAMD) of, for example, novel refrigerants or in drug discovery, where their ability to accurately predict physical properties for new species aids in evaluating the performance of a hypothetical molecule. Here, the assignment of groups must be done thousands of times and, in some cases, for rather complex molecules. This is the primary motivator for the development of GCIdentifier.

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