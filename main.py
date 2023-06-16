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


def get_element_number(atom_coordinates, num_columns):
    x = atom_coordinates[0]
    y = atom_coordinates[1]
    z = atom_coordinates[2]
    return z * num_columns ** 2 + y * num_columns + x + 1


def convert_xyz_to_mesh(filepath, show_initial, show_min_max_coord, number_of_nodes_in_each_line):
    molecule = read(filepath)
    print(f" [{molecule.get_chemical_formula()}]", "--> ", end="")
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

    new_positions = positions - min_coordinates
    # print(f"new_positions after the movement to make everything in the positive quadrant: \n{new_positions}")
    # print(f"cube length:{actual_cube_length}")

    cube_length_scaling_factor = number_of_nodes_in_each_line / actual_cube_length
    # print(f"scaling_factor:{cube_length_scaling_factor}")
    new_cube_length = actual_cube_length * cube_length_scaling_factor
    # print(f"new cube length:{new_cube_length}")
    new_positions = new_positions * cube_length_scaling_factor
    # print(f"final position of atoms after the scaling: \n{new_positions}")
    rounded_positions_to_match_with_mesh_grids = new_positions.astype(int)
    # print(f"final position of atoms after rounding: \n{rounded_positions_to_match_with_mesh_grids}")
    u, c = np.unique(rounded_positions_to_match_with_mesh_grids, return_counts=True)
    assert len(u[c > 1]) != 0
    molecule.set_positions(rounded_positions_to_match_with_mesh_grids)
    result = []
    for atom in molecule:
        result.append(f"{atom.symbol}{int(get_element_number(atom.position, number_of_nodes_in_each_line))}")
    result.append(f"-{round(actual_cube_length, 6)}-{number_of_nodes_in_each_line}")
    result_str = "".join(result)
    print(result_str)
    return result_str


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
    before_is_number = False
    for i in range(len(mesh_part)):
        c = mesh_part[i]
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
                temp.append(c)
            else:
                numbers.append("".join(temp))
                temp.clear()
                temp.append(c)
            before_is_number = False
    print(symbols)
    print(numbers)
    print(len(symbols), len(numbers))


def convert_mesh_to_xyz(mesh_str: str):
    mesh_part, cube_length, mesh_number_in_line = mesh_str.split("-")
    print(mesh_part)
    print(cube_length)
    print(mesh_number_in_line)
    convert_mesh_str_to_sym_and_num(mesh_part)


xyz_files = os.listdir("./xyz_fragments")
for index, xyz_filename in enumerate(xyz_files):
    xyz_filepath = f"./xyz_fragments/{xyz_filename}"
    print(xyz_filename.replace(".xyz", ""), end="")
    mesh_str = convert_xyz_to_mesh(xyz_filepath, show_initial=False, show_min_max_coord=False, number_of_nodes_in_each_line=100)
    convert_mesh_to_xyz(mesh_str)
    exit()
