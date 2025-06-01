'''
This example shows how to visualize natural orbitals using PySCF's cubeview module.'''

from pyscf import gto, scf, mcscf, cubeview

if __name__=="__main__":
    mol = gto.M(
    atom = """
 C     0.      0.      0.658 
 C     0.      0.     -0.658 
 H     0.      0.924   1.252 
 H     0.     -0.924   1.252 
 H     0.     -0.924  -1.252 
 H     0.      0.924  -1.252 
""",
    basis = "6-31g*",
    spin = 0, charge = 0, symmetry = "D2h",
    verbose = 3)

    mf = scf.RHF(mol)
    mf.xc = 'b3lyp'
    mf.kernel()

    mc =  mcscf.CASSCF(mf, 2, 2)
    mc.natorb = True
    mc.kernel()

    cv = cubeview.viewer(mc)
    cv.prepare(mo_list=[7,8,9,10]) 
    cv.serve()
