````@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: GCIdentifier.jl
  text: Automated molecular fragmentation for group-contribution methods
  image:
    src: /logo.png
    alt: GCIdentifier
  tagline: Assign any SMILES to UNIFAC, Joback, SAFT-γ Mie, gcPC-SAFT and more — or propose entirely new groups
  actions:
    - theme: brand
      text: Getting started
      link: /group_search
    - theme: alt
      text: View on GitHub
      link: https://github.com/ClapeyronThermo/GCIdentifier.jl

features:
  - icon: 🧩
    title: Group Assignment
    details: Fragment any molecule into groups for ten supported group-contribution methods
    link: /group_search

  - icon: 🔍
    title: Find Missing Groups
    details: Automatically propose new groups for molecules not yet covered by existing methods
    link: /missing_groups

  - icon: ⚗️
    title: Custom Groups
    details: Define your own SMARTS-based groups and plug them into any supported method
    link: /custom_groups
---
````

```@meta
CurrentModule = GCIdentifier
```

GCIdentifier.jl fragments a molecular SMILES (or name) into the groups defined by existing group-contribution methods — UNIFAC, Joback, SAFT-γ Mie, gcPC-SAFT and others — and can automatically propose new groups for molecules not yet covered. It is designed for high-throughput use in computer-aided molecular design (CAMD).

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
