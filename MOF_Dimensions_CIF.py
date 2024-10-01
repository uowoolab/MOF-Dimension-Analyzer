#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2020-01-21: created by S.B.Wiggin, the Cambridge Crystallographic Data Centre
#
# 2024-09-30: Modified by Olivier Marchand, University of Ottawa, to support CIF structures as input.

"""
    MOF_Dimensions_CIF.py  - processes CIF structures by performing two expansions of a polymeric 
    network and produces minimum area bounding boxes. Comparison of the two bounding boxes gives the growth dimensions 
    of the framework.
"""

from ccdc.io import EntryReader
import numpy as np
import argparse
import csv
import os.path
import time


def generate_bounding_box(atoms):

    all_pts = np.array([[c for c in atom.coordinates] for atom in atoms])
    cov = np.cov(all_pts, rowvar=False)
    evals, evecs = np.linalg.eig(cov)

    lengths = np.sqrt(evals)
    lengths += 0.0001  # Adjust for exceedingly small box dimensions
    lengths.sort()

    return lengths


def dimensionality(entry):
    """
    Calculates the dimensionality of the crystal.
    Grows 2 instances of the polymeric unit (four cycles and seven cycles) and calculates the change in size
    The number of cycles needs to grow a representative part of the framework, but more cycles will take longer
    """

    t0 = time.time()
    shell1 = entry.crystal.polymer_expansion(repetitions=4)
    shell2 = entry.crystal.polymer_expansion(repetitions=7)

    t1 = time.time()
    time_taken = str(t1 - t0)
    print("total time elapsed for polymer expansion is % s" % time_taken)

    lengths1 = generate_bounding_box(shell1.atoms)
    print("number of atoms in first expansion " + str(len(shell1.atoms)))
    s1, m1, l1 = lengths1

    lengths2 = generate_bounding_box(shell2.atoms)
    print("number of atoms in second expansion " + str(len(shell2.atoms)))
    s2, m2, l2 = lengths2

    ratios = [s2 / s1, m2 / m1, l2 / l1]
    print("ratio of box dimensions from first and second expansions: " + str(ratios))
    thresh = 1.15
    ndims = [r > thresh for r in ratios].count(True)
    return ndims


def analyse_structures(user_gcd_input, writer, n_structures):
    """
    Analyzes structures from a CIF file and appends the results to a CSV file.

    Parameters:
    - user_gcd_input: The path to the CIF file to be analyzed.
    - writer: A csv.writer object for writing the analysis results to a CSV file.
    - n_structures: The current count of structures processed, used for indexing.

    Returns:
    - n_structures: Updated count of structures processed.
    - n_mof: Count of MOF structures identified in the current CIF file.
    - n_non_mof: Count of non-MOF structures identified in the current CIF file.
    """

    n_mof = 0
    n_non_mof = 0

    try:
        csd_reader = EntryReader(user_gcd_input, format="cif")

        for entry in csd_reader:
            print("\nCIF file: " + str(entry.identifier))
            n_structures += 1
            count_polymers = 0
            for component in entry.molecule.components:
                if component.is_polymeric:
                    count_polymers += 1
            if count_polymers > 1:
                print("Multiple polymer units present")
            if entry.molecule.heaviest_component.is_polymeric:
                n_mof += 1
                framework = entry.molecule.heaviest_component
                framework.remove_hydrogens()

                entry.crystal.molecule = framework

                fig = dimensionality(entry)

                if fig == 0:
                    dimension = "0D non-MOF"
                elif fig == 1:
                    dimension = "1D chain"
                elif fig == 2:
                    dimension = "2D sheet"
                elif fig == 3:
                    dimension = "3D framework"
            else:
                n_non_mof += 1
                dimension = "No polymeric bonds detected"

            print(
                "Framework dimensions for CIF file %s: %s"
                % (entry.identifier, dimension)
            )
            writer.writerow((entry.identifier, dimension, n_structures))

    except Exception as e:
        file = os.path.split(user_gcd_input)[1]
        error_message = str(e).strip()
        print(f"\nFailed to process CIF file {file}:\n{error_message}")
        writer.writerow((file, "FAILED", n_structures))
        n_structures += 1

    return n_structures, n_mof, n_non_mof


def get_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", help="Absolute path of your database")
    parser.add_argument("-o", "--output", help="Results CSV filename")

    return parser.parse_args()


def main():

    args = get_args()

    n_structures = 0
    if args.output is None or args.input is None:
        return

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("Filename", "Dimensionality", "File Order"))

        for cif_file in os.listdir(args.input):
            cif_path = os.path.join(args.input, cif_file)
            if os.path.isfile(cif_path):
                n_structures, n_mof, n_non_mof = analyse_structures(
                    cif_path, writer, n_structures
                )


if __name__ == "__main__":
    main()