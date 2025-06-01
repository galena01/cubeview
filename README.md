CubeView - PySCF-based Orbital Visualization Toolkit with 3Dmol.js Integration
=========================

2025-05-31

* Version 0.1

![](https://raw.githubusercontent.com/galena01/cubeview/refs/heads/main/webpage.jpg)

Generates cube files using PySCF's cubegen module and visualizes molecular orbitals through 3Dmol.js. Key displayed properties include:

- Orbital symmetry
- Orbital energy
- Occupation numbers

The Python script produces fully static HTML files that can be:
* Exported to directories
* Viewed in local browsers
* Embedded in web platforms (e.g., github.io pages)

The static nature of the webpage also reduces security risks from potential attacks.

Currently supported orbital types:

* Canonical MO from RHF(RKS)/UHF(UKS)
* Natural orbitals from CASSCF
* Localized MOs

Installation
-------
Install to Python site-packages:

```
pip install git+https://github.com/galena01/cubeview
```

Usage
-------

Visualize Molecular Orbitals

```python
from pyscf import cubeview

# Perform calculation ...

cv = cubeview.viewer(mf)
cv.prepare(mo_list=[1,2,5,6],ao_component=1)    # May consume significant storage space when visualizing numerous orbitals

cv.serve()  # To start a local web server

cv.export('./xxx') # To export the visualization to a directory
```

Other examples can be found in the examples directory.
