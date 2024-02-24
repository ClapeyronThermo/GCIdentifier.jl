---
title: 'GCIdentifier.jl: A Julia package for identifying molecular fragments from SMILES'
tags:
  - Julia
  - Group-contribution
  - Thermodynamics
  - Molecular Design
authors:
  - name: Pierre J. Walker
    orcid: 0000-0001-8628-6561
    corresponding: true
    affiliation: "1, 2"
  - name: Andrés Riedemann
    corresponding: true
    affiliation: 3
  - name: Zhen-Gang Wang
    affiliation: 1
affiliations:
 - name: Division of Chemistry and Chemical Engineering, California Institute of Technology, Pasadena, California 91125, United States
   index: 1
 - name: Department of Chemical Engineering, Imperial College, London SW7 2AZ, United Kingdom
   index: 2
 - name: Departamento de Ingeniería Química, Universidad de Concepción, Concepción 4030000, Chile
   index: 3
date: 9 February 2024
bibliography: paper.bib
---

# Summary
GCIdentifier.jl is an open-source toolkit for the automatic identification of group fragments based on the name of a molecule or its SMILES. Obtaining chemical properties of species, such as heat capacities [@jobackEstimationPureComponentProperties1987] or solvation free energies [@plattsEstimationMolecularLinear2000], will typically involve a set of parameters which represent a given species. Unfortunately, in these cases, the parameters obtained only apply to a specific species and cannot be transferred to others. An alternative approach would be to split the molecule into moieties, known as groups, each of which will have their own parameters associated with them. The combination of these groups (and their associated parameters) will represent the entire molecule. The benefit is that these groups can be combined in different ways such that they represent a new molecule. This type of approach is known as group contribution, where multiple examples of such approaches exist [@weidlichModifiedUNIFACModel1987;@walkerNewPredictiveGroupContribution2020;@chungGroupContributionMachine2022;@papaioannouGroupContributionMethodology2014]. An example of this process is shown in figure 1.

![Fragmentation of ibuprofen into UNIFAC groups.](figures/ibuprofen.pdf)

Unfortunately, the challenge with using group-contribution approaches is the assignment of the groups to represent a given species. While this assignment can be done manually, it is more-convenient and, as discussed later, efficient to automate this process. Indeed, this is the exact objective of GCIdentifier. By simply feeding a species name or SMILES, along with the group-contribution approach one wishes to use, the group assignment is done automatically:
```julia
using GCIdentifier, ChemicalIdentifiers

groups = get_groups_from_name("ibuprofen", UNIFACGroups)
```
The output from this function can then be used in other packages, such as Clapeyron [@walkerClapeyronJlExtensible2022], to obtain chemical properties.


# Statement of need
Group-contribution approaches are vital when it comes to computer-aided molecular design (CAMD) of, for example, novel refrigerants [@sahinidisDesignAlternativeRefrigerants2003] or in drug discovery[@houADMEEvaluationDrug2004]. Here, the assignment of groups must be done thousands of times and, in some cases, for rather complex molecules. This is the primary motivator for the development of GCIdentifier. While other packages with similar functionalities have been developed in other languages, GCIdentifier stands apart for multiple reasons.

GCIdentifier is the first of such packages to be compatible with multiple group-contribution approaches, such as UNIFAC and SAFT-$\gamma$ Mie. By standardising the representation of groups using SMARTS and leveraging the powerful MolecularGraph [@seiji_matsuoka_2024_10478701] package, our group-identification code can be used with any existing group-contribution thermodynamic model. This extends to group-contribution approaches which require information about the connectivity between groups [@sauerComparisonHomoHeterosegmented2014] where, by simply specifying `connectivity=true` within the `get_groups_from_name` function, the connectivity matrix between groups will automatically be generated.

While packages in other languages are able to generate groups from _existing_ group databases, GCIdentifiers is able to systematically propose _new_ groups for a given molecule. Consider a case where an existing group-contribution framework is unable to cover all atoms present in a molecule. GCIdentifier is able to consider these un-represented atoms and propose a list of new groups. From this list, users will be able to determine which groups they should obtain new parameters for. In the extreme case where we wish to generate a list of all possible groups that represent a molecule, GCIdentifier will automatically split the molecule into groups, from which either the user or a set of built-in heuristics can then decide which set best represent the molecule. 

These two features present within GCIdentifier have potential applications beyond thermodynamic modelling, such as the development of molecular dynamics forcefields which could be integrated into packages such as Molly [@greenerDifferentiableSimulationDevelop2023].

# Acknowledgments
Z-G.W. acknowledges funding from Hong Kong Quantum AI Lab, AIR\@InnoHK of the Hong Kong Government.

# References