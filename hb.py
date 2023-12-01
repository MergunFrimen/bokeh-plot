import math
import numpy as np
from Bio import PDB
from bokeh.plotting import figure, show
from bokeh.models import ColorBar
from bokeh.transform import linear_cmap
from bokeh.palettes import Viridis256, gray
from bokeh.models import ColumnDataSource
from urllib.request import urlretrieve


def calculate_distance_matrix(structure):
    model = structure[0]
    num_residues = len(list(model.get_residues()))
    distance_matrix = np.zeros((num_residues, num_residues))

    for chain in model:
        for res1 in chain.get_residues():
            for res2 in chain.get_residues():
                ca1 = get_alpha_carbon(res1)
                ca2 = get_alpha_carbon(res2)
                if ca1 is not None and ca2 is not None:
                    res1_index = res1.id[1] - 1
                    res2_index = res2.id[1] - 1
                    distance = np.linalg.norm(ca1 - ca2)
                    distance_matrix[res1_index, res2_index] = distance
                    distance_matrix[res2_index, res1_index] = distance

    return num_residues, distance_matrix


def get_alpha_carbon(residue):
    for atom in residue.get_atoms():
        if atom.id == 'CA':
            return atom.get_vector()
    return None


def plot_contact_map(num_residues, distance_matrix, threshold=8.0):
    p = figure(
                width=800,
                height=800,
                x_range=(1, num_residues),
                y_range=(1, num_residues),
                x_axis_label="Amino acid index",
                y_axis_label="Amino acid index",
                sizing_mode='scale_width',
                x_axis_location="below",
                y_axis_location="left",
                toolbar_location=None,
                tools="hover"
               )

    x, y = np.indices(distance_matrix.shape)
    x = x.flatten()
    y = y.flatten()
    distances = distance_matrix.flatten()

    source = ColumnDataSource(data=dict(x=x, y=y, value=distances))
    mapper = linear_cmap(field_name='value', palette=gray(2), low=0, high=threshold)

    p.square(
        'x',
        'y',
        source=source,
        size=15,
        line_color=None,
        fill_alpha=1,
        color=mapper
        )

    p.hover.tooltips = [('Residue', '@y, @x'), ('Distance', '@value')]

    show(p)


if __name__ == "__main__":
    url = 'https://files.rcsb.org/download/1AH9.pdb'
    pdb_filename = "1ah9.pdb"
    urlretrieve(url, pdb_filename)
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_filename)

    num_residues, distance_matrix = calculate_distance_matrix(structure)

    plot_contact_map(num_residues, distance_matrix, 15.0)
