#
# BLAntarctic v0.1 - 'muscle' module
#
# Hace MUSCLE, instancia los resultados. Genera árbol filogenético y plotea.
#

import os
from subprocess import PIPE
from subprocess import Popen
import matplotlib
import matplotlib.pyplot as plt

from Bio.Blast import NCBIXML
from Bio.Align.Applications import MuscleCommandline
from Bio import Phylo
from Bio import SeqIO

def multi_align(estruct_dir):
    """
    Hace un alineamiento múltiple de las secuencias contenidas en cada
    multifasta de la carpeta creada por 'blastp.build_multifasta_for_muscle'.
    Almacena los resultados en la carpeta 'aligned_matches'.
    """

    # Crea el directorio de salida del alineamiento.
    os.mkdir(estruct_dir.file_in_results_dir("aligned_matches"))

    # Recorre los fichero del directoiro 'unaligned_matches'.
    for filename in os.listdir(estruct_dir.file_in_results_dir("/unaligned_matches")):

        try:
            base_name = "{}_{}".format(filename.split("_")[0], filename.split("_")[1])

            muscle_cline = MuscleCommandline(
                    input = estruct_dir.file_in_results_dir("/unaligned_matches/{}".format(filename)),
                    out = estruct_dir.file_in_results_dir("/aligned_matches/{}_aligned.fa".format(base_name)))

            stdout, stderr = muscle_cline()

        except:
            # .DS_Store, etc.
            pass

class Query:
    """
    Sus instancias almacenan información sobre los queries utilizados, necesaria
    para graficar más tarde.
    """

    def __init__(self, id, aligned_seq):
        self.id = id
        self.aligned_seq = aligned_seq
        self.non_gapped = []
        self.aligned_domains = []

    def ilustrate_alignment(self, start, end):
        """
        Recoge las regiones de la secuencia que no son gaps en el alineamiento
        (que no son guiones).
        """
        self.non_gapped.append((start,end))

    def add_aligned_domain(self, accession, start, end):
        """
        Recoge las posiciones de los dominios conservados en relación a la
        secuencia alineada (incluyendo guiones).
        """
        self.aligned_domains.append((accession, start, end))


def parse_alignment(estruct_dir, proteins):
    """
    Parsea el resultado del alineamiento, e incorpora a las instancias de la
    clase Protein las secuencias con los gaps necesarios para cuadrar el
    alineamiento.
    """

    path = estruct_dir.file_in_results_dir("aligned_matches")

    queries = []

    for filename in os.listdir(path):
        # Recorre los archivos del directorio.
        query_id = filename.split("_")[0] + "_" + filename.split("_")[1]
        with open(path+"/"+filename) as input_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                for protein in proteins:
                    if protein.query_id == query_id:
                        if protein.subject_id == record.id:
                            protein.add_alignment(record.seq)

                if record.id == query_id:
                    queries.append(Query(record.id, record.seq))
    return(queries)

def build_tree(estruct_dir):
    """
    Construye árbol filogenético en formato NW.
    """

    # Crea el directorio de salida.
    os.mkdir(estruct_dir.file_in_results_dir("trees_nw"))

    # Recorre los fichero del directoiro 'unaligned_matches'.
    for filename in os.listdir(estruct_dir.file_in_results_dir("/aligned_matches")):

        try:

            input = estruct_dir.file_in_results_dir("/aligned_matches/{}".format(filename))

            base_name = filename.split("_")[0]

            output = estruct_dir.file_in_results_dir("/trees_nw/{}_tree.nw".format(base_name))

            process = Popen(
                    ['muscle', '-in', input, '-out', output, '-maketree', '-cluster', 'neighborjoining'],
                    stdout=PIPE,
                    stderr=PIPE)

            stdout, stderr = process.communicate()

        except:
            # .DS_Store, etc.
            pass

def plot_tree(estruct_dir):
    """
    Genera una representación gráfica del árbol en formato PDF a partir del
    fichero .nw que crea MUSCLE.
    """

    source = estruct_dir.file_in_results_dir("trees_nw")
    os.mkdir(estruct_dir.file_in_results_dir("trees_plot"))

    for filename in os.listdir(source):

        out = estruct_dir.file_in_results_dir("trees_plot")
        base_name = filename.split("_")[0]

        # Genera el gráfico con Phylo.
        tree = Phylo.read("{}/{}".format(source, filename), 'newick')
        tree.ladderize()

        matplotlib.rc('font', size=6)

        fig = plt.figure(figsize=(10, 20), dpi=100)

        # Inserta el árbol en un gráfico de Matplotlib.
        axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=axes, do_show=False)
        plt.title("{} BLAST matches phylogeny".format(base_name))

        # Lo guarda en PDF.
        plt.savefig("{}/{}_tree_plot.pdf".format(out, base_name), dpi=100)

        plt.close()
