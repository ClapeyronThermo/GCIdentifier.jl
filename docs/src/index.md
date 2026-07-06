````@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: GCIdentifier.jl
  text: Automated molecular fragmentation for group-contribution methods
  image:
    src: logo.png
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

## Citing GCIdentifier.jl

If you are using GCIdentifier for your research work, please cite the following:

```text
@article{GCIdentifier-2024,
	title = {GCIdentifier.jl: A Julia package for identifying molecular fragments from SMILES},
	author = {Walker, Pierre J. and Riedemann, Andrés and Wang, Zhen-Gang},
	journal = {J. Open Source Softw.},
	volume = {9},
  number = {96},
	year = {2024},
	pages = {6453},
	url = {https://joss.theoj.org/papers/10.21105/joss.06453},
	doi = {10.21105/joss.06453},
}
```

## Related packages

````@raw html
<div class="related-pkg-grid">
  <a class="related-pkg-card" href="https://clapeyronthermo.github.io/Clapeyron.jl/" target="_blank" rel="noreferrer">
    <div class="related-pkg-logo-wrap">
      <img class="related-pkg-logo" src="/assets/related_clapeyron_logo.svg" alt="Clapeyron.jl" />
    </div>
    <h3 class="related-pkg-title">Clapeyron.jl</h3>
    <p class="related-pkg-details">Provides every bulk equation of state cDFT builds its inhomogeneous functionals on top of, and is required alongside cDFT for essentially all use.</p>
  </a>
  <a class="related-pkg-card" href="https://clapeyronthermo.github.io/Langmuir.jl/dev/" target="_blank" rel="noreferrer">
    <div class="related-pkg-logo-wrap">
      <img class="related-pkg-logo" src="/assets/related_langmuir_logo.png" alt="Langmuir.jl" />
    </div>
    <h3 class="related-pkg-title">Langmuir.jl</h3>
    <p class="related-pkg-details">Single- and multi-component adsorption equilibrium models, complementary to cDFT's own adsorption isotherm calculations.</p>
  </a>
</div>

<style>
.related-pkg-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
  gap: 16px;
  margin: 16px 0 32px;
}

.related-pkg-card {
  display: flex;
  flex-direction: column;
  align-items: center;
  text-align: center;
  padding: 24px;
  border: 1px solid var(--vp-c-bg-soft);
  border-radius: 12px;
  background-color: var(--vp-c-bg-soft);
  text-decoration: none !important;
  transition: border-color 0.25s, background-color 0.25s;
}

.related-pkg-card:hover {
  border-color: var(--vp-c-brand-1);
}

.related-pkg-logo-wrap {
  display: flex;
  align-items: center;
  justify-content: center;
  width: 100%;
  height: 96px;
  margin-bottom: 8px;
}

.related-pkg-logo {
  max-height: 100%;
  max-width: 100%;
  width: auto;
  height: auto;
  object-fit: contain;
}

.related-pkg-title {
  margin: 0;
  line-height: 24px;
  font-size: 16px;
  font-weight: 600;
  color: var(--vp-c-text-1);
  border-top: none;
  padding-top: 0;
}

.related-pkg-details {
  flex-grow: 1;
  margin: 8px 0 0;
  line-height: 22px;
  font-size: 14px;
  font-weight: 500;
  color: var(--vp-c-text-2);
}
</style>
````