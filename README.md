# Berea_PNM

Primary drainage in static pore network modelling: a comparative study

Research output: Thesis › Master's Thesis

Authors

Omidreza Amrollahinasab Mahdiabad

Chair of Reservoir Engineering
Abstract

Darcy's law or other similar approaches are implemented to model the flow at the continuum scale. However, the accuracy of these models depends on porous media parameters such as porosity, relative permeability, capillary pressure relationships, etc (Bhattad et al., 2011). One option is to obtain these parameters experimentally, but in recent years, the development of rock imaging tools that image the porous rock and the fluids inside them has improved our understanding of the flow in porous media. In combination with publicly available numerical tools, it enables us to simulate the flow in these media (J. Blunt, 2017). Pore network modeling (PNM) is a widely used technique to simulate multiphase transport in porous materials, and one of the open-source available software packages to do so is OpenPNM. Several porous media research groups have developed it as a pore network modeling tool. In pore network modeling, a pore network is extracted from reconstructed microtomography images and the corresponding numerical calculations are done to predict the transport properties of the porous media (Gostick J. A., 2016). Pore network modeling leads to the calculation of the petrophysical and advanced rock flow properties. In this work, the OpenPNM python package is implemented and used to simulate the primary drainage capillary pressure and relative permeability from microtomography images of sandstone rocks. Pore networks are extracted using the Porespy python package which uses an algorithm named SNOW to extract pore network from segmented microtomography rock images. The model is validated using the experimental results of capillary pressure, absolute permeability, and primary drainage relative permeabilities from well-known rock samples. Simulation results have proved to match the experimentally obtained results within an acceptable range of uncertainty. A sensitivity analysis is then conducted to understand the influence of modeling parameters, interfacial tension, and contact angle on the results. Finally, a workflow is presented to predict the petrophysical and transport properties of unknown rock samples. The implementation of the workflow is then demonstrated on a sandstone rock of interest.

ref: https://pure.unileoben.ac.at/portal/en/publications/primary-drainage-in-static-pore-network-modelling-a-comparative-study(bb925aa8-102d-4239-878c-e644c074cb8f).html?customType=theses

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/omidreza-amrollahi/Berea_PNM/HEAD)
