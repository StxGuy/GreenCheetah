# GreenCheetah

![Logo](https://github.com/StxGuy/GreenCheetah/blob/main/gsheetah.svg)

## About

    GreenCheetah: A Non-Equilibrium Green's Function (NEGF) approach for quantum transport in Fortran.

[![Generic badge](https://img.shields.io/badge/GitHub-StxGuy/GreenCheetah-<COLOR>.svg)](https://github.com/StxGuy/GreenCheetah)

## Usage

```use GreenCheetah```

[...]

```
type(GF) :: G
   
G = GF(Potential,delta_y,delta_x)
```   
   
##### Transmission

``` 
t = G%Transmission(number_of_points,maximum_energy)
```

##### Density of states

``` 
D = G%DOS(energy)
``` 


## Credits


    @misc{daCunha2021,
        author       = {C. R. da Cunha},
        title        = {{GreenCheetah, a library for non-equilibrium Green's functions}},
        month        = apr,
        year         = 2021,
        version      = {1.0},
        publisher    = {GitHub},
        url          = {https://github.com/StxGuy/GreenCheetah}
        }
        
## License

Copyright (c) 2021 - Carlo R. da Cunha, "GreenCheetah"
<carlo.requiao@gmail.com>
