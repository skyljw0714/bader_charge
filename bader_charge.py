import sys
import numpy as np
import ase
import ase.io
import os
import ase.neighborlist
import networkx as nx 


class BCharge:
    def __init__(self):
        self.get_files() 
        self.atomic_charge = self.get_atomic_charge() 
        self.charge_by_atom = self.get_charge_by_atom()
    
    def _get_chargesum(self):
        os.system('chgsum.pl AECCAR0 AECCAR2')

    def _get_acf(self):
        if not os.path.exists('./CHGCAR_sum'):
            self._get_chargesum()
        os.system('bader CHGCAR -ref CHGCAR_sum')
        
    def get_files(self):
        folder = '.'
        if not os.path.exists(folder + '/ACF.dat'):
            self._get_acf()  # make acf file only if it does not exist
        self.pos = folder + '/POSCAR'
        self.pot = folder + '/POTCAR'
        self.acf = folder + '/ACF.dat' 

    #make graph 
    def _get_covalent_neighbor_list(self, atoms, SCALE=1.2):
        cutoffs = ase.neighborlist.natural_cutoffs(atoms)
        cutoffs = [SCALE * c for c in cutoffs]
        return ase.neighborlist.neighbor_list("ijD", atoms, cutoff=cutoffs)

    #make graph 
    def get_graph(self, atoms):
        i, j, _ = self._get_covalent_neighbor_list(atoms)
        G = nx.Graph()
        G.add_edges_from(zip(i, j))
        return G

    def get_connected_atoms(self):
        atoms = ase.io.read(self.pos)
        graph = self.get_graph(atoms) 
        subgraphs = nx.connected_components(graph) 
        components = []
        symbols = [] 
        
        for i, c in enumerate(list(subgraphs)):
            index = list(c)
            symbol = sorted(list(set(atoms[index].get_chemical_symbols())))
    
            components.append(index)
            symbols.append(symbol)
            
            if symbol == ['N', 'O']: #NO2 
                assert len(index) == 3     
        assert ['Mo', 'S'] in symbols
        return components, symbols

    def get_connected_atoms_charge(self):
        charge = self.atomic_charge
        components, symbols = self.get_connected_atoms()
        for c,s in zip(components, symbols):    
            print(c)
            print(charge[c])
            print(sum(charge[c]))
            print(s)

    def get_atoms_from_poscar(self): #read POSCAR 
        with open(self.pos, 'r') as f:
            data = f.readlines() 
        atoms = get_items_from_list(data[5].split(' '))
        atoms = remove_n_from_list(atoms)
        nums = list(map(int, get_items_from_list(data[6].split(' '))))
        return {a : [n] for a,n in zip(atoms,nums)} #{'Mo': 9, 'S': 18}

    def get_zval_from_potcar(self):
        atoms = self.get_atoms_from_poscar() 
        with open(self.pot, 'r') as f:
            data = f.readlines()
        for i,d in enumerate(data):
            if '_PBE' in d and 'TITEL' not in d: #start line 
                atom = d.split(' ')[3].split('_')[0]
                zval = float(data[i+1])
                #print(d)
                #print(zval)
                atoms[atom].append(zval)
        return atoms
    
    def get_charge_from_acf(self):
        with open(self.acf, 'r') as f:
            data = f.readlines()
        charge = [] 
        for d in data[2:-4]:
            charge.append(float(get_items_from_list(d.split(' '))[4])) #charge from acf.dat
        return np.array(charge)
        
    def get_atomic_charge(self):
        atoms = self.get_zval_from_potcar() 
        charge = self.get_charge_from_acf() 
        zval = []
        for n, c in atoms.values():
            zval += [c] * n 
        zval = np.array(zval) 
        assert len(charge) == len(zval)
        return zval - charge
    
    def get_charge_by_atom(self):
        charge_by_atom = {}
        i = 0
        for atom, number in self.get_atoms_from_poscar().items():
            charge_by_atom[atom] = self.atomic_charge[i:i+number[0]]
            i = i+number[0]
        return charge_by_atom

def get_items_from_list(lst): #remove '' items in lst
    new = []
    for l in lst:
        if (l == '') or (l == '\n'):
            continue
        new.append(l)
    return new

def remove_n_from_list(lst): #remove \n in lst items 
    return [i.replace('\n','') for i in lst]  

if __name__ == "__main__":
    bc = BCharge()
    bc.get_connected_atoms_charge() 
    #charge = BCharge().get_atomic_charge()
    #print(charge)

