import numpy as np
import ovito
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys
from dataclasses import dataclass
from itertools import product, combinations_with_replacement


@dataclass
class SROModifier:

    type_map: dict = None
    
    def __call__(self, frame: int, data: ovito.data.DataCollection):
    
        # get bond topology, i.e. a list of (a, b) pairs
        # bonds are not double counted, need to be careful
    
        topology = np.array(data.particles.bonds.topology[...], dtype=int)
        
        # get types of each particle, and assosciated quantities
        
        types = np.array(data.particles['Particle Type'][...], dtype=int)
        num_atoms = len(types)
        unique_types = list(set(types))
        num_types = len(unique_types)
        
        # calculate concentrations
        
        concentrations = np.zeros(num_types)
        for i, type_ in enumerate(unique_types):
            concentrations[i] = np.mean(types == type_)
        
        # initialize an array that stores the number of (i, j) bonds we have
        # where i and j are particle types
        
        bond_types_count_array = np.zeros((num_types, num_types))
        
        for index, bond in enumerate(topology):
        
            # get particle indexes from bond
        
            first_atom, second_atom = bond
            
            # get types of those particle indexes
            
            first_type, second_type = types[first_atom], types[second_atom]
            
            # count that bond
            
            bond_types_count_array[first_type - 1, second_type - 1] += 1
            
        # normalize by the number of bonds we have, make array symmetric
            
        probability_array = bond_types_count_array / len(topology)
        probability_array += probability_array.T
        
        # calculate sro_array
        # outer gives array such that O[i, j] = c[i] * c[j]
        
        sro_array = 1.0 - probability_array / np.outer(concentrations, concentrations)
        
        # create iterable of type pairs
        
        iterable = product(unique_types, repeat=2)
        for i, j in iterable:
        
            # if a type map is provided, store values with atom names
            # else, store with integer type labels
        
            sro = sro_array[i - 1, j - 1]
            if not self.type_map:
                key = f'sro_{i}{j}'
            else:
                key = f'sro_{self.type_map[i]}{self.type_map[j]}'
            data.attributes[key] = sro
            
            
def tuple_to_pair(t: tuple):

    t1, t2 = t
    
    return f'{t1}{t2}'
            
        
def main():
    
    input_dump_file = sys.argv[1]
    output_plot_file = sys.argv[2]
    
    # define k-nn cutoffs, i.e. a bond is in k-nn shell if first_cutoff < bond_length < second_cutoff
    
    first_cutoff = 0.0
    second_cutoff = 3.2
    
    # create ovito pipeline, add a bonds modifier to create needed topology

    pipeline = ovito.io.import_file(input_dump_file)
    bonds_modifier = ovito.modifiers.CreateBondsModifier(lower_cutoff=first_cutoff, cutoff=second_cutoff)
    pipeline.modifiers.append(bonds_modifier)
    
    # define a type map, mapping integer labels to atom types
    
    type_map = {1: 'Fe', 2: 'Ni', 3: 'Cr', 4: 'Co', 5: 'Mn'}
    
    # create SRO modifier which calculates all SRO's at each timestep
    
    sro_modifier = SROModifier(type_map=type_map)
    pipeline.modifiers.append(sro_modifier)
    
    # initialize atom pairs
    pairs = list(combinations_with_replacement(type_map.values(), 2))
    
    # initialize arrays to store values in
    
    frames = np.arange(pipeline.source.num_frames)
    timestep = np.zeros(frames.shape)
    sro_params = np.zeros((pipeline.source.num_frames, len(pairs)))
    
    for frame in frames:
    
        data = pipeline.compute(frame)
        timestep[frame] = data.attributes['Timestep']
        for pair_index, pair in enumerate(pairs):
            sro_params[frame, pair_index] = data.attributes[f'sro_{tuple_to_pair(pair)}']
        
    # 1e-6 assumes metal units and a 1 fs timestep, femtoseconds to nanoseconds
    
    for pair_index, pair in enumerate(pairs):
    
        plt.scatter(1e-6 * timestep, sro_params[:, pair_index], zorder=6, label=tuple_to_pair(pair), edgecolor='black', alpha=0.5)
    plt.legend(title='$ij$')
    plt.grid()
    plt.xlabel('time (ns)')
    plt.ylabel('first $ij$ Cowley SRO parameter')
    plt.savefig(output_plot_file, dpi=800, bbox_inches='tight')
    
    
if __name__ == '__main__':

    main()

