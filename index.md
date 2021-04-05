## Welcome to GreenCheetah Documentation

![Logo](https://github.com/StxGuy/GreenCheetah/blob/main/gsheetah.svg)

GreenCheetah is a Fortran library for a Non-Equilibrium Green's function (NEGF) approach for solving problems of quantum transport.


## Installation

make

## Usage

To initialize:

    use GreenCheetah

    [...]
    
    type(GreeFunction) :: G

    G = GreenFunction(Potential,number_of_rows,number_of_columns)

#### Transmission Function
    t = G%Transmission(number_of_points,maximum_energy)

#### Density of States
    DOS = G%DOS(energy)

## Credits
    Projet GreenCheetah https://github.com/StxGuy/GreenCheetah
    Copyright (c) 2021 Carlo R. da Cunha (carlo.requiao@gmail.com)

    @misc{daCunha2021,
      author = {C. R. da Cunha},
      title = {GreenCheetah},
      year = {2021},
      publisher = {GitHub},
      version = {1.0},
      journal = {GitHub repository},
      url= {https://github.com/StxGuy/GreenCheetah}
    }


