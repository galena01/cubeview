import tempfile
from pyscf.tools import cubegen
from pyscf import symm, scf, mcscf, gto
from pyscf.gto.mole import tostring
import os, sys, subprocess, shutil
import json
import gzip 

def process_ao(coeff, labels, n_ao):
    contributions = [(abs(c)**2, idx) for idx, c in enumerate(coeff)]
    contributions.sort(reverse=True)
    top_2 = contributions[0:n_ao]
    ao_info = ", ".join([f"{labels[idx]} ({contrib:.4f})" for contrib, idx in top_2])
    return ao_info

def gen_mo_list_json(mol, mo_coeff, mo_energy, mo_occ, irreps, ao_component=0, spin=''):
    '''
    Generate a JSON object of Molecular Orbitals with their energies and symmetry labels.
    Parameters
    mol : Mole
        Molecule object.
    mo_coeff : ndarray
        MO coefficients.
    ao_component : int, optional
        Number of AO components to show. Default is 0.
    spin: str, optional
        Spin component ('a','b','') to show. Default is ''.
    '''

    n_orb = mo_coeff.shape[0]
    if mo_energy is None:
        mo_energy = [0.0] * n_orb
    if mo_occ is None:
        mo_occ = [0.0] * n_orb
    if irreps is None:
        irreps = ["A"] * n_orb

    assert(ao_component>=0) # Number of components must be positive
    assert(ao_component<=mo_coeff.shape[1]) # Number of components must be less than or equal to the number of MOs.
    assert(mo_coeff.shape[0]==mol.nao_nr()) # MO coefficients must have the same number of rows as the number of AOs.


    irrep_counts = {}
    
    data = []

    for i in range(1,n_orb+1):
        energy = mo_energy[i-1]
        occ = mo_occ[i-1]

        irrep_name = irreps[i-1]
        irrep_counts[irrep_name] = irrep_counts.get(irrep_name, 0) + 1
        symm_text = f"{irrep_name}.{irrep_counts[irrep_name]}"

        # idx = f'{i}/{spin}' if spin else f'{i}'
        
        entry = {
            "index": i,
            "irrep": symm_text,
            "spin": spin,
            "occ": round(occ, 4),
            "energy": round(energy, 4),
        }
        
        if ao_component>0:
            ao_info = process_ao(mo_coeff[:,i-1], mol.ao_labels(), ao_component)
            entry["ao_components"] = ao_info
        
        data.append(entry)

    return data

def gen_cube_file(mol, mo_coeff, fn, resolution=0.3, margin=5.0, ftmp="/dev/shm/tmp.cube"):
    """Common method to generate cube files and compress them"""
    cubegen.orbital(mol, ftmp, mo_coeff, resolution=resolution, margin=margin)
    # Compress cube file and delete original file
    with open(ftmp, 'rb') as f_in:
        with gzip.open(fn, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(ftmp)

class CubeViewBase:
    def __init__(self, mf_or_mol, **kwargs):
        self._workdir = tempfile.TemporaryDirectory()
        self.workdir = self._workdir.name
        print(f"Temporary directory created at {self.workdir}")
        # self.mf = mf
        if isinstance(mf_or_mol, gto.MoleBase):
            self.mol = mf_or_mol
        else:
            self.mf = mf_or_mol
        
    def prepare(self, **kwargs):
        '''Common preparation steps for all viewers'''
        template_src = os.path.join(os.path.dirname(__file__), 'static')
        # Copy all files from the template directory to the workdir
        if not os.path.exists(self.workdir):
            raise FileNotFoundError(f"Work directory {self.workdir} does not exist.")
        
        for item in os.listdir(template_src):
            s = os.path.join(template_src, item)
            d = os.path.join(self.workdir, item)
            if os.path.isdir(s):
                shutil.copytree(s, d, dirs_exist_ok=True)
            else:
                shutil.copy2(s, d)
    def serve(self, port=8000, **kwargs):
        '''Start HTTP server in temp dir to serve interactive cube viewer'''


        python_executable = sys.executable
        subprocess.run([python_executable, "-m", "http.server", str(port)], cwd=self.workdir, **kwargs)


    def cleanup(self):
        '''Cleanup temporary files'''
        self._workdir.cleanup()

    def export(self, path):
        '''Export the generated files to a specified path'''
        if not os.path.exists(path):
            os.makedirs(path)
        # recursively copy the workdir to the specified path
        for item in os.listdir(self.workdir):
            s = os.path.join(self.workdir, item)
            d = os.path.join(path, item)
            if os.path.isdir(s):
                shutil.copytree(s, d, dirs_exist_ok=True)
            else:
                shutil.copy2(s, d)

        # Create a README file
        with open(os.path.join(path, "README.txt"), "w") as f:
            f.write("This directory contains the generated cube files and metadata for molecular orbitals.\n")
            f.write("You can visualize the orbitals using the provided HTML viewer.\n")
            f.write("To start the viewer, run `python -m http.server [port]` in this directory and open your browser to http://localhost:[port]\n")


class RHFCubeView(CubeViewBase):
    def prepare(self, mo_list=None, **kwargs):
        super().prepare(**kwargs)
        
        irreps = symm.label_orb_symm(self.mf.mol, self.mf.mol.irrep_name, self.mf.mol.symm_orb, self.mf.mo_coeff) if self.mf.mol.symmetry else None
        # Generate JSON file
        with open(f"{self.workdir}/orbitals.json", "w") as f:
            json.dump(
                gen_mo_list_json(
                    self.mf.mol, self.mf.mo_coeff, self.mf.mo_energy, self.mf.mo_occ,
                    irreps, 
                    kwargs.get("ao_component", 0)),
                f, indent=4
            )

        # Generate xyz file
        open(f"{self.workdir}/mol.xyz",'w').write(
            tostring(self.mf.mol,format='xyz')
        )

        # Generate cube files
        if mo_list is None:
            mo_list = range(1, self.mf.mo_coeff.shape[1]+1)
        os.makedirs(f"{self.workdir}/cubes", exist_ok=True)
        for i in mo_list:
            print(f"Generating cube file for MO {i}")
            fn = f"{self.workdir}/cubes/{i}.cube.gz"
            gen_cube_file(self.mf.mol, self.mf.mo_coeff[:,i-1], fn)

class UHFCubeView(CubeViewBase):

    def prepare(self, mo_list=None, **kwargs):
        super().prepare(**kwargs)
        
        irreps_a = symm.label_orb_symm(self.mf.mol, self.mf.mol.irrep_name, self.mf.mol.symm_orb, self.mf.mo_coeff[0]) if self.mf.mol.symmetry else None
        irreps_b = symm.label_orb_symm(self.mf.mol, self.mf.mol.irrep_name, self.mf.mol.symm_orb, self.mf.mo_coeff[1]) if self.mf.mol.symmetry else None
        # Generate JSON file
        with open(f"{self.workdir}/orbitals.json", "w") as f:
            json.dump(
                gen_mo_list_json(
                    self.mf.mol, self.mf.mo_coeff[0], self.mf.mo_energy[0], self.mf.mo_occ[0], 
                    irreps_a,
                    kwargs.get("ao_component", 0), spin='a') + \
                gen_mo_list_json(
                    self.mf.mol, self.mf.mo_coeff[1], self.mf.mo_energy[1], self.mf.mo_occ[1],
                    irreps_b,
                    kwargs.get("ao_component", 0), spin='b'),
                f, indent=4
            )

        # Generate xyz file
        open(f"{self.workdir}/mol.xyz",'w').write(
            tostring(self.mf.mol,format='xyz')
        )

        n_orb = self.mf.mo_coeff[0].shape[1]
        if mo_list is None:
            mo_list = [range(1, n_orb+1), range(1, n_orb+1)]
        os.makedirs(f"{self.workdir}/cubes", exist_ok=True)
        for spin_id, spin in enumerate(['alpha', 'beta']):
            for i in mo_list[spin_id]:
                print(f"Generating cube file for {spin} MO {i}")
                fn = f"{self.workdir}/cubes/{spin}_{i}.cube.gz"
                gen_cube_file(self.mf.mol, self.mf.mo_coeff[spin_id][:,i-1], fn)

class CASSCFCubeView(RHFCubeView):
    pass

class CustomCubeView(CubeViewBase):
    
    def __init__(self, mol, mo_coeff=None, mo_energy=None, mo_occ=None, irreps=None):
        super().__init__(mol)
        
        self.mo_coeff = mo_coeff
        self.mo_energy = mo_energy
        self.mo_occ = mo_occ
        self.irreps = irreps
        assert(mo_coeff is not None), "MO coefficients must be provided for CustomCubeView."

    def prepare(self, mo_list=None, **kwargs):
        super().prepare(**kwargs)
        
        # Generate JSON file
        with open(f"{self.workdir}/orbitals.json", "w") as f:
            json.dump(
                gen_mo_list_json(self.mol, self.mo_coeff, self.mo_energy, self.mo_occ,
                                 self.irreps, kwargs.get("ao_component", 0)),
                f, indent=4
            )

        # Generate xyz file
        open(f"{self.workdir}/mol.xyz",'w').write(
            tostring(self.mol,format='xyz')
        )

        # Generate cube files
        if mo_list is None:
            mo_list = range(1, self.mo_coeff.shape[1]+1)
        os.makedirs(f"{self.workdir}/cubes", exist_ok=True)
        for i in mo_list:
            print(f"Generating cube file for MO {i}")
            fn = f"{self.workdir}/cubes/{i}.cube.gz"
            gen_cube_file(self.mol, self.mo_coeff[:,i-1], fn)

def viewer(mf_or_mol, mo_coeff=None, mo_energy=None, mo_occ=None,irreps=None):
    '''cube view factory'''
    if isinstance(mf_or_mol, gto.MoleBase):
        return CustomCubeView(mf_or_mol, mo_coeff=mo_coeff, mo_energy=mo_energy, mo_occ=mo_occ, irreps=irreps)
    
    if isinstance(mf_or_mol, scf.uhf.UHF):
        return UHFCubeView(mf_or_mol)
    elif isinstance(mf_or_mol, scf.rhf.RHF):
        return RHFCubeView(mf_or_mol)
    elif isinstance(mf_or_mol, mcscf.mc1step.CASSCF):   # casscf
        return CASSCFCubeView(mf_or_mol)
    else:
        raise ValueError(f"Unsupported method type: {type(mf_or_mol)}.")


