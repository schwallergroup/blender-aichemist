## 🤓 Course overview and learning outcomes 

The goal of this mini-workshop is to give you a brief introduction on how to use python to generate atomistic structures ready for Blender. We’ll also provide you with materials for further learning and a few ideas to get you started. 🚀

## Before we continue...

To run this notebook/code snipets you will need only a few packages installed (if you don't have them already). In additon, we will use autoadsorbate (wich only requires rdkit and ase):
```python
pip install rdkit
pip install ase
pip install autoadsorbate
```
We will use these basic packages to construct structures of our choice, using ```SMILES``` and ```*SMILES``` (for more info see: https://github.com/basf/autoadsorbate). We will start by making a few molecular structures, we will construct complex geometries, and finally we will make a trajectory so that we can render images in bulk.

We will use free and portable (no install required) software:
Blender (download: )
Inkscape (download: )

### Init
Let us start with importing the required packages:

```python
from ase import Atoms
from ase.visualize.plot import plot_atoms
from ase.io import read, write
import autoadsorbate as au
```

#### Converting SMILES into 3D

We can easily convert any ```SMILES``` string into a molecular geometry using the following code:

```python
f = au.Fragment('Cc1ccccc1', to_initialize=1)
atoms = f.get_conformer(0)
plot_atoms(atoms, rotation='-60x')
```    
![png](getting_started_files/getting_started_3_1.png)


If we have multiple molecules we can use a simple loop to collect all molecules in a list:

```python
smiles = [
    'Cc1ccccc1',
    'Brc1ccccc1',
    '[Mg](C)Br',
    '[Mg](Br)Br',
]

molecules = []
for s in smiles:
    molecules.append(au.Fragment(s, to_initialize=1).get_conformer(0)) 
```

### Metal Complex structures

Unfortunately, rdkit cannot help us with the geometries of the metal complex structures. But it is relatively easy to prepare structures for illustration purposes (or as input for e.g. quantum chemical simulations) using a few tricks:

#### Helper Function

AutoAdsorbate (https://github.com/basf/autoadsorbate) comes with a convenient function ```align_to_vector``` that will orient a ligand (provided as ```*SMILES*```), to an arbitrary cartesian vector. We only need to provide the central metal atom, surrogate smiles of the species attachec to it, and the vectors to define the complex geometry.

```python
def get_complex_structure(central_atom, smiles, vectors):
    """"
    helper function that creates a complex (/interediate) out of
    a central metal atom, complex geometry and *SMILES (surrogate smiles).
    For more information on *SMILES see: https://github.com/basf/autoadsorbate
    """
    complex = central_atom.copy()
    
    for i, v in enumerate(vectors):
        # prepare ligand
        ligand = au.Fragment(
            smiles[i], to_initialize=1
            ).get_conformer(0)                        # converts *SMILES string to XYZ oriented towards Z
        ligand = au.Smile.align_to_vector(ligand, v)  # aligns the conformer to the gemetry of the complex
        del ligand[0]                                 # removes the surrogate atom
        complex+=ligand
    return complex
```

#### make ligand with *SMILES


```python
# '*[P](C)(C)C' we use Cl as surrogate atom: 'Cl[P](C)(C)C'

f = au.Fragment('Cl[P](C)(C)C', to_initialize=1) #here we are using the Cl-P bond to orient the ligand
atoms = f.get_conformer(0)
plot_atoms(atoms, rotation='-60x')
```
    
![png](getting_started_files/getting_started_9_2.png)



```python
# https://en.wikipedia.org/wiki/Tetrahedron
complex_geometry = {
    'linear': [[0,0,1], [0,0,-1]],
    'tetrahedron': [[1,1,1],[-1,-1,1],[1,-1,-1],[-1,1,-1]]
}
```

#### PdL2


```python
Pd_atom = Atoms(symbols=['Pd'], positions=[[0,0,0]])

smiles = [
    'Cl[P](C)(C)C',
    'Cl[P](C)(C)C',
]

vectors = complex_geometry['linear']

PdL2 = get_complex_structure(Pd_atom, smiles, vectors)

plot_atoms(PdL2, rotation='-60x')
```

    [20:13:32] UFFTYPER: Unrecognized charge state for atom: 1
    [20:13:32] UFFTYPER: Unrecognized charge state for atom: 1





    <Axes: >




    
![png](getting_started_files/getting_started_12_2.png)
    


#### PdL2Y2


```python
Pd_atom = Atoms(symbols=['Pd'], positions=[[0,0,0]])

smiles = [
    ['Cl[P](C)(C)C','Cl[P](C)(C)C', 'ClBr', 'Clc1ccccc1'],
    ['Cl[P](C)(C)C','Cl[P](C)(C)C', 'ClC', 'Clc1ccccc1']
]

vectors = complex_geometry['tetrahedron']

Pd_comps = []

for s in smiles:
    Pd_comps.append(get_complex_structure(Pd_atom, s, vectors))

view(Pd_comps)
```

    [19:59:32] UFFTYPER: Unrecognized charge state for atom: 1
    [19:59:32] UFFTYPER: Unrecognized charge state for atom: 1
    [19:59:32] UFFTYPER: Unrecognized charge state for atom: 1
    [19:59:32] UFFTYPER: Unrecognized charge state for atom: 1





    <Popen: returncode: None args: ['/home/djrm/venv/mace_env/bin/python', '-m',...>




```python
collected_traj = molecules + [PdL2] + Pd_comps
write('./collected_traj.xyz', collected_traj)
```

### Make blenderized traj


```python
collected_traj = read('./collected_traj.xyz', index=':')

for a in collected_traj:
    print(a.get_chemical_formula(empirical=True))
```

    C7H8
    C6H5Br
    CH3BrMg
    Br2Mg
    C6H18P2Pd
    C12H23BrP2Pd
    C13H26P2Pd



```python
blend_traj = au.utils.get_blenderized(collected_traj, hide_spot=[0,0,-100])

for a in blend_traj:
    print(a.get_chemical_formula(empirical=True))
    
write('./collected_traj.xyz', collected_traj)
```

    C13H26Br2MgP2Pd
    C13H26Br2MgP2Pd
    C13H26Br2MgP2Pd
    C13H26Br2MgP2Pd
    C13H26Br2MgP2Pd
    C13H26Br2MgP2Pd
    C13H26Br2MgP2Pd



```python

```
