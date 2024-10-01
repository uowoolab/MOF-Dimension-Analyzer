#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2020-01-21: created by S.B.Wiggin, the Cambridge Crystallographic Data Centre
#
# 2024-09-30: Modified by Olivier Marchand, University of Ottawa, to support CIF structures as input
# and parallelized execution to enhance efficiency and speed.

"""
    MOF_Dimensions_CIF_Parallel.py  - processes CIF structures by performing two expansions of a polymeric 
    network and produces minimum area bounding boxes. Comparison of the two bounding boxes gives the growth dimensions 
    of the framework. The code now supports parallel execution for increased efficiency and speed.
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
from ccdc.io import EntryReader
import numpy as np
import argparse
import csv
import os.path
import time


class ResultCollector:
    """
    A class to collect information as the code runs and print the results at the end.
    """

    def __init__(self, identifier):
        self.identifier = identifier
        self.results = []

    def add_info(self, info):
        """Adds information to the results."""
        self.results.append(info)

    def print_results(self):
        """Prints all accumulated results."""
        print(f"\nCIF file: {self.identifier}")
        for result in self.results:
            print(result)


def generate_bounding_box(atoms):

    all_pts = np.array([[c for c in atom.coordinates] for atom in atoms])
    cov = np.cov(all_pts, rowvar=False)
    evals, evecs = np.linalg.eig(cov)

    lengths = np.sqrt(evals)
    lengths += 0.0001  # Adjust for exceedingly small box dimensions
    lengths.sort()

    return lengths


def dimensionality(entry, result_collector):
    """
    Calculates the dimensionality of the crystal.
    Grows 2 instances of the polymeric unit (four cycles and seven cycles) and calculates the change in size.
    """

    t0 = time.time()
    shell1 = entry.crystal.polymer_expansion(repetitions=4)
    shell2 = entry.crystal.polymer_expansion(repetitions=7)

    t1 = time.time()
    time_taken = str(t1 - t0)
    result_collector.add_info(
        f"Total time elapsed for polymer expansion: {time_taken} seconds"
    )

    lengths1 = generate_bounding_box(shell1.atoms)
    result_collector.add_info(
        f"Number of atoms in first expansion: {len(shell1.atoms)}"
    )

    lengths2 = generate_bounding_box(shell2.atoms)
    result_collector.add_info(
        f"Number of atoms in second expansion: {len(shell2.atoms)}"
    )

    ratios = [s2 / s1 for s1, s2 in zip(lengths1, lengths2)]
    result_collector.add_info(
        f"Ratio of box dimensions from first and second expansions: {ratios}"
    )

    thresh = 1.15
    ndims = sum(r > thresh for r in ratios)

    return ndims


def analyse_structure(cif_path):
    """
    Analyzes structures from a CIF file and appends the results to a CSV file.

    Parameters:
    - cif_path: The path to the CIF file to be analyzed.

    Returns:
    - results: A list of tuples containing the CIF file identifier and the dimension.
    """

    n_mof, n_non_mof = 0, 0
    results = []

    try:
        csd_reader = EntryReader(cif_path, format="cif")

        for entry in csd_reader:
            dimension = "No polymeric bonds detected"
            result_collector = ResultCollector(entry.identifier)

            if any(component.is_polymeric for component in entry.molecule.components):
                n_mof += 1
                framework = entry.molecule.heaviest_component
                framework.remove_hydrogens()

                entry.crystal.molecule = framework

                ndims = dimensionality(entry, result_collector)

                if ndims == 0:
                    dimension = "0D non-MOF"
                elif ndims == 1:
                    dimension = "1D chain"
                elif ndims == 2:
                    dimension = "2D sheet"
                elif ndims == 3:
                    dimension = "3D framework"
                result_collector.add_info(f"Framework dimension: {dimension}")

            else:
                n_non_mof += 1

            result_collector.print_results()
            results.append((entry.identifier, dimension))

    except Exception as e:
        error_message = str(e).strip()
        print(
            f"\nFailed to process CIF file {os.path.basename(cif_path)}:\n{error_message}"
        )
        results.append((os.path.basename(cif_path), "FAILED"))

    return results


def get_args():
    parser = argparse.ArgumentParser(description="Parallel MOF Dimensionality Analysis")
    parser.add_argument(
        "-i", "--input", required=True, help="Directory containing CIF files"
    )
    parser.add_argument("-o", "--output", required=True, help="Output CSV file path")
    parser.add_argument(
        "-n",
        "--num_cpus",
        type=int,
        default=1,
        help="Number of CPUs to use for processing. Defaults to 1 if unspecified.",
    )
    return parser.parse_args()


def main():
    args = get_args()

    cif_files = [
        os.path.join(args.input, f)
        for f in os.listdir(args.input)
        if os.path.isfile(os.path.join(args.input, f))
    ]

    with open(args.output, "w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["Filename", "Dimensionality"])

        with ProcessPoolExecutor(max_workers=args.num_cpus) as executor:
            futures = {
                executor.submit(analyse_structure, cif_path): cif_path
                for cif_path in cif_files
            }

            for future in as_completed(futures):
                results = future.result()
                for result in results:
                    writer.writerow(result)
                csv_file.flush()


if __name__ == "__main__":
    main()