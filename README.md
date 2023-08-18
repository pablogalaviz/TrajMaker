# TrajMaker

Takes an input structure in POSCAR format and creates an XDATCAR trajectory with a prescribed oscillatory dynamics.


## Getting Started

git clone repository 
```shell
git clone https://github.com/pablogalaviz/TrajMaker.git
cd TrajMaker
```

### Prerequisites

This project needs Python 3. Create a new venv:
```shell
python3 -m venv venv
```
Activate the virtual environment:
```shell
source venv/bin/activate
```

## Installation
Installation can be done with the following command:
```shell
pip install .
```
## Execution
The script `trajmaker.py` needs a configuration file in yaml format.
```shell
trajmaker.py -c input.yaml
```
See for example `test/input.yaml`.

### Units
The atoms move around their initial position by a displacement given by:

```latex
dᵢ = Sᵢ(I,J,K) Aᵢ sin (2 π (ω t - ϕᵢ))
```
For `i ∈ {x,y,z}`. Here `I,J,K` denotes the supercell index, `S` is a scale factor for each supercell, `A` 
the amplitude of the oscillation in angstrom. 
The frequency `ω` and the time `t` are both defined in the configuration file. 
Therefore, the units are reciprocal, if `final_time` is considered in picoseconds, `ω` is in THz.



## History

First release Aug 2023

## Credits

Author(s): 


Pablo Galaviz | galavizp@ansto.gov.au


## License

TrajMaker is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

TrajMaker is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with TrajMaker.  

See <http://www.gnu.org/licenses/>.
