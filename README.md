# cif2hkl
A program that computes structure factors |F^2| for neutrons, x-rays, and electrons

**cif2hkl** - CrysFML based utility for generating .lau files for McStas.

Action: Read a CIF/CFL/SHX/PCR crystallographic description and generates a HKL F^2 reflection list.

Compilation
-----------
```
gfortran -O2 -o cif2hkl cif2hkl.F90 -lm 
```

Syntax
------
```bash
cif2hkl --xtal file.cif cif2hkl --powder file.cif
```

Use
---

The cif2hkl program is a command line tool, that uses arguments. The syntax and help are obtained with:
```bash
cif2hkl --help
```
It reads a file with crystal structure information, and computes the |F^2| structure factors. The supported input file formats are:

| Extension | Format | Link |
|-----------|--------|------|
| CIF | Crystallographic Information File | https://en.wikipedia.org/wiki/Crystallographic_Information_File |
| PCR | FullProf control input file       | https://www.ill.eu/sites/fullprof/ |
| CFL | CrysFML input file                | https://code.ill.fr/scientific-software/crysfml |
| SHX | ShelX input file                  | https://shelx.uni-goettingen.de/ |
| INS | ShelX instruction file            | https://shelx.uni-goettingen.de/ |
| RES | ShelX result                      | https://shelx.uni-goettingen.de/ |

The general syntax is:
```
cif2hkl [options][-o outfile] file1 file2 ...
```
and can treat files in series. The result is a file with a readable header and reflection list with columns
```
[ H K L Multiplicity Sin(Theta/Lambda) d_spacing |F|^2 ]
```
which is used by e.g. [McStas](http://www.mcstas.org) neutron ray-tracing software. However, as X-rays and electronic scattering structure factors can also be computed, the tool has a wider use.

The available options on the comand line are:
```
--help     or -h    Show this help
--version  or -v    Display program version
--out FILE          Specify the name of the next output file.
   -o FILE            Default is to add .hkl to the initial file name.
--lambda LAMBDA     Set the incoming probe wavelength [Angs].
   -l    LAMBDA       Default is 0.5
--powder   or -p    Generate a list of unique HKL reflections (for powders). Default.
--xtal     or -x    Generate a list of all HKL reflections (for single crystals).
--verbose           Display processing details.
--no-outout-files   Just read the CIF/CFL/ShellX file (for checking).
```

Example: cif2hkl -o CaF2.laz CaF2.cfl

Credits and License
-------
This software is (c) E. Farhi, Synchrotron Soleil with license EUPL (GPL2 compatible).
Part of the iFit <http://ifit.mccode.org> suite.
It is based on CrysFML (CFML) available at <https://code.ill.fr/scientific-software/crysfml>, but all required modules are all included in the cif2hkl source code.
