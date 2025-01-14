## Issues

## Functionality that needs to be completed

GUI:
- Getting Started
- About

## Enhancements

- Database of known atomic data (displacement energy, FP energy, density of elements)

- GUI: In material def form, add a button to calc density based on composition (by simple mixing)

- Provide progress info for HDF5 i/o operations

- Make user-defined tallies for various events. E.g.
  - Implantation (position, atomic species)
  - Vacancy (position, atomic species)
  - Ion escape (have to distinguish backscattered/transmitted ions)
  - Recoil (PKA or other)
    - energy / damage energy
    - position
    - atom
    - vacancies generated
  - The tallies have their own mesh which can be defined in rect, spherical, cylindrical coordinates 

- Handle surface effects