This solver is built on top of the standard `pimpleFoam` framework and extends it by incorporating a passive scalar transport equation. It can be used to trace how a scalar field mixes under the influence of turbulence (passive).
- Includes adaptive mesh refinement (AMR) guided by a user-defined `gradU` field, which can be modified to suit specific needs.
