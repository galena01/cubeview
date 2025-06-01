'''
This example shows how to export molecular orbitals to a directory.
'''

from pyscf import gto, scf, cubeview

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
    cv = cubeview.viewer(mf)
    cv.prepare(mo_list=[1,2,5,6])
    cv.export('./mycube/')
