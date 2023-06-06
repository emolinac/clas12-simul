Modification of LEPTO 6.5.1
===========================
 
Monte Carlo Event Generator.

It depends on `pythia5721.f` and `jetset7410.f`. The original program is called `lepto-6.5.1.f`. The modification is described in `qp1.f`.

## Requirements

* [**ROOT**](https://root.cern/)

* [**cernlib**](http://github.com/utfsm-eg2-data-analysis/cernlib)

* [**clas_software**](http://github.com/utfsm-eg2-data-analysis/clas_software)

## Compilation

Compile by running `make`.

## Usage

To execute the compiled binary file `lepto.exe` (located in `bin/`), there must exist a text file called `lepto.txt` in the same directory with the following content:

```
<Nevts> <tarA> <tarZ> <pid>
```

This means that `<Nevts>` events will be generated with at least one particle `<pid>`, simulating interactions with a nuclear target with A = `<tarA>` and Z = `<tarZ>`.

For example, for a quick test of 12 events containing omega particles from DIS with Deuterium target, one could execute:

```
echo "12 2 1 223" > lepto.txt
./lepto.exe
```

## References

1. T. Sjöstrand. “High-energy-physics event generation with Pythia 5.7 and Jetset 7.4” (1994). DOI:[10.1016/0010-4655(94)90132-5](http://doi.org/10.1016/0010-4655(94)90132-5).

2. G. Ingelman *et al*. “Lepto 6.5 - A Monte Carlo generator for deep inelastic lepton-nucleon scattering” (1997). DOI:[10.1016/S0010-4655(96)00157-9](http://doi.org/10.1016/S0010-4655(96)00157-9).
