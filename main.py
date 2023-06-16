import os
import shutil
import xml.etree.ElementTree as ET
from ase import Atoms
from ase.io import read
from ase.optimize import BFGS
from ase.calculators.emt import EMT
import numpy as np
from bs4 import BeautifulSoup
from time import sleep
from tabulate import tabulate

def get_element_number(atom_coordinates, num_columns):
    i = atom_coordinates[0]
    j = atom_coordinates[1]
    k = atom_coordinates[2]
    return i * num_columns ** 2 + j * num_columns + k + 1


def convert_xyz_to_mesh(filepath, show_initial, show_min_max_coord, number_of_nodes_in_each_line):
    molecule = read(filepath)
    print(molecule.get_chemical_formula())
    if show_initial:
        print(f"initial positions: \n{molecule.get_positions()}")
    positions = molecule.get_positions()
    T_positions = positions.T
    max_coordinates = np.array([np.max(T_positions[i]) for i in range(3)])
    min_coordinates = np.array([np.min(T_positions[i]) for i in range(3)])
    if show_min_max_coord:
        print(f"minimum point of the cube coordinates:{min_coordinates}")
        print(f"maximum point of the cube coordinates:{max_coordinates}")
    actual_cube_length = np.max(max_coordinates - min_coordinates)

    new_positions1 = positions - min_coordinates
    # print(f"new_positions after the movement to make everything in the positive quadrant: \n{tabulate(new_positions1, headers=['X', 'Y', 'Z'])}")
    # print(f"cube length:{actual_cube_length}")

    cube_length_scaling_factor = (number_of_nodes_in_each_line - 1) / actual_cube_length
    # print(f"scaling_factor:{cube_length_scaling_factor}")
    new_cube_length = actual_cube_length * cube_length_scaling_factor
    # print(f"new cube length:{new_cube_length}")
    new_positions2 = new_positions1 * cube_length_scaling_factor
    # print(f"final position of atoms after the scaling: \n{new_positions}")
    rounded_positions_to_match_with_mesh_grids = new_positions2.astype(int)
    # print(tabulate(rounded_positions_to_match_with_mesh_grids, headers=['i', 'j', 'k']))
    # print(f"final position of atoms after rounding: \n{rounded_positions_to_match_with_mesh_grids}")
    u, c = np.unique(rounded_positions_to_match_with_mesh_grids, return_counts=True)
    assert len(u[c > 1]) != 0
    molecule.set_positions(rounded_positions_to_match_with_mesh_grids)
    result = []
    for atom in molecule:
        result.append(f"{atom.symbol}{int(get_element_number(atom.position, number_of_nodes_in_each_line))}")
    result.append(f"-{round(actual_cube_length, 6)}-{number_of_nodes_in_each_line}")
    result_str = "".join(result)
    return result_str, new_positions1


# convert_xyz_to_mesh("water.xyz", show_initial=True, show_min_max_coord=True, number_of_nodes_in_each_line=100)
# convert_xyz_to_mesh("water.xyz", show_initial=True, show_min_max_coord=True, number_of_nodes_in_each_line=1000)


def get_all_cml_files(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            print(file_path)
            src = file_path
            des = "C:\\Users\\Amin\\PycharmProjects\\pythonProject\\cml_fragments\\{}".format(file)
            # shutil.copyfile(file_path, des)


def convert_cml_to_xyz(cml_file, xyz_file):
    # Read the CML file
    with open(cml_file, 'r') as f:
        cml_data = f.read()

    # Parse the CML file with BeautifulSoup
    soup = BeautifulSoup(cml_data, 'lxml')

    # Extract atomic coordinates and elements
    atoms = soup.find_all('atom')
    coordinates = []
    elements = []
    for atom in atoms:
        coordinates.append([float(atom['x3']), float(atom['y3']), float(atom['z3'])])
        elements.append(atom['elementtype'])
    # Write the XYZ file
    with open(xyz_file, 'w') as f:
        # Write the number of atoms and a comment line
        f.write(str(len(coordinates)) + '\n')
        f.write('Converted from CML to XYZ\n')

        # Write atomic coordinates with element symbols
        for i in range(len(coordinates)):
            line = f"{elements[i]} {coordinates[i][0]} {coordinates[i][1]} {coordinates[i][2]}\n"
            f.write(line)


def convert_mesh_str_to_sym_and_num(mesh_part: str):
    symbols = []
    numbers = []
    temp = []
    before_is_number = None
    for i in range(len(mesh_part)):
        c = mesh_part[i]
        if before_is_number is None:
            temp.append(c)
            before_is_number = False
            continue
        if c.isdigit():
            if before_is_number:
                temp.append(c)
            else:
                #     the previous one was symbol
                symbols.append("".join(temp))
                temp.clear()
                temp.append(c)

            before_is_number = True
        else:
            if before_is_number:
                numbers.append("".join(temp))
                temp.clear()
                temp.append(c)

            else:
                temp.append(c)
            before_is_number = False
    numbers.append("".join(temp))
    temp.clear()
    return symbols, numbers

def get_indices(element_number, number_nodes_in_line):
    element_number -= 1
    k = element_number % number_nodes_in_line
    j = (element_number // number_nodes_in_line) % number_nodes_in_line
    i = element_number // (number_nodes_in_line * number_nodes_in_line)
    return i, j, k




def convert_mesh_to_xyz(mesh_str: str):
    mesh_part, cube_length, mesh_number_in_line = mesh_str.split("-")
    symbols, numbers = convert_mesh_str_to_sym_and_num(mesh_part)
    ijks = [get_indices(int(number), int(mesh_number_in_line)) for number in numbers]
    distance_between_two_nodes = float(cube_length) / (int(mesh_number_in_line) - 1)
    xyzs = [(ijk[0] * distance_between_two_nodes, ijk[1] * distance_between_two_nodes, ijk[2] * distance_between_two_nodes) for ijk in ijks]
    return ijks, xyzs
    # return symbols, numbers, cube_length, mesh_number_in_line


xyz_files = os.listdir("./xyz_fragments")
for index, xyz_filename in enumerate(xyz_files):
    xyz_filepath = f"./xyz_fragments/{xyz_filename}"
    print(xyz_filename.replace(".xyz", ""))
    mesh_str, initial_positions = convert_xyz_to_mesh(xyz_filepath, show_initial=False, show_min_max_coord=False, number_of_nodes_in_each_line=1000)
    print(mesh_str)
    ijks_from_mesh, xyzs_from_mesh = convert_mesh_to_xyz(mesh_str)
    result_table_headers = ["i", "j", "k", "actual X", "actual Y", "actual Z", "converted X", "converted Y", "converted Z", "% error**2"]
    # print(ijks_from_mesh)
    # print(initial_positions)
    # print(xyzs_from_mesh)
    result_table = [(ijk_from_mesh[0], ijk_from_mesh[1], ijk_from_mesh[2], act_xyz[0], act_xyz[1], act_xyz[2], xyz[0], xyz[1], xyz[2], (np.sqrt(((xyz[0] - act_xyz[0])**2 + (xyz[1] - act_xyz[1])**2 + (xyz[2] - act_xyz[2])**2)/(act_xyz[0]**2 + act_xyz[1]**2 + act_xyz[2]**2)) * 100)) for ijk_from_mesh, act_xyz, xyz in zip(ijks_from_mesh, initial_positions, xyzs_from_mesh)]
    print(tabulate(result_table, headers=result_table_headers))
    print(100*"-")
