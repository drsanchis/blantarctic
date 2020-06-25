#
# BLAntarctic v0.1 - 'blastp' module
#
# Hace BLAST, filtra los resultados, los almacena en instancias, genera ficheros
# output.
#

from Bio.Blast.Applications import NcbiblastpCommandline
import os

from Bio import Seq
from Bio import SeqIO

def blast_it(estruct_dir, evalue, coverage):
    """blastp

    Sobre el conjunto de subjects que se indique, se hace un blastp de la
    secuencia query. Se genera un fichero de salida con un nombre contruido a
    partir del argumento 'nombre_proyecto'.
    """

    subject = estruct_dir.file_in_data_dir('genomes_multifasta.fa')
    query = estruct_dir.file_in_data_dir('queries.fa')
    out = estruct_dir.file_in_results_dir("blast_results.tsv")

    blastp_cline = NcbiblastpCommandline(
            query = query,
            subject = subject,
            outfmt = '"6 qseqid sseqid sseq evalue qcovs pident qstart qend"',
            out = out,
            evalue = evalue,
            qcov_hsp_perc = coverage)

    stdout, stderr = blastp_cline()

def filter_ident(estruct_dir, ident):

    input = estruct_dir.file_in_results_dir("blast_results.tsv")
    output = estruct_dir.file_in_results_dir("filtered_blast_results.tsv")

    with open(output, 'a') as out:
        with open(input, 'r') as file:
            for line in file:
                split_line = line.split("\t")
                if float(split_line[5]) > ident:
                    out.write(line)

class Protein:
    queries = []

    def __init__(self, idx,
                 query_id,subj_id, subj_seq,
                 eval, query_cov, ident, qstart, qend):
        """
        Inicializa la instancia. Almacena en variables de instancia la
        información correspondiente a cada resultado del BLAST.
        """
        self.index = int(idx) # Asigna ID numerico a la instancia.
        self.query_id = query_id
        self.subject_id = subj_id
        self.subject_seq = subj_seq
        self.evalue = float(eval)
        self.query_cov = float(query_cov)
        self.ident = float(ident)
        self.range = (qstart, qend)

        # Separa el subject id en locus_tag y especie.
        split_id = subj_id.split("@")
        self.subject_accs = str(split_id[0])
        self.subject_spec = str(split_id[1])
        self.subject_id_spaced = "{} @ {}".format(
                self.subject_accs,
                self.subject_spec)

        # Inicializa estas listas, sobre las cuales se harán appends después.
        self.domains = []
        self.aligned_domains = []
        self.non_gapped = []

        # Lleva un registro de todos los queries para los cuales se han
        # guardado resultados.
        Protein.queries.append(query_id)
        Protein.temp1 = set(Protein.queries)
        Protein.queries = list(Protein.temp1)

    def add_domain(self, accession, start, end):
        """
        Añade un dominio conservado detectado.
        """
        self.domains.append((accession, start, end))

    def add_alignment(self, aligned_seq):
        """
        Añade a la instancia la secuencia correspondiente a esa proteína dentro
        del alineamiento múltiple (con gaps) para después graficarlo.
        """
        self.aligned_seq = aligned_seq

    def add_aligned_domain(self, accession, start, end):
        """
        Añade un dominio conservado, pero tomando como referencia las posiciones
        dentro de la secuencia alineada (con gaps), en lugar de en la secuencia
        original, para después graficar.
        """
        self.aligned_domains.append((accession, start, end))

    def ilustrate_alignment(self, start, end):
        """
        Añade regiones de la secuenia alineada que se corresponden con
        aminoácidos y no con gaps (guiones), para después graficar.
        """
        self.non_gapped.append((start,end))

    def import_queries_list(queries_list):
        """
        Permite recuperar la lista de queries al abrir el fichero *.bapj.
        """
        Protein.queries = queries_list


def parse_blast_results(estruct_dir, Protein_class):
    """
    Parsea el fichero tsv con los resultados del BLAST, y crea una instancia de
    la clase 'Protein' para cada match, incluyendo en los atributos de
    dicha instancia la información relevante del match.
    """
    proteins = [] # Genera una lista que contendrá todas las instancias.

    with open (estruct_dir.file_in_results_dir("filtered_blast_results.tsv"), 'r') as results:
        for idx, line in enumerate(results, start=1): # Recorre líneas del tsv
            campos = line.split("\t")

            qseqid = campos[1]

            # Recupera la secuencia original del match, no la secuencia
            # recortada que aparece en los resultados del BLAST.
            with open(estruct_dir.file_in_data_dir("genomes_multifasta.fa")) as filehandle:
                for record in SeqIO.parse(filehandle, "fasta"):
                    if qseqid == record.id:
                        original_sequence = str(record.seq)

            proteins.append(Protein_class(idx, campos[0],
                                          campos[1], original_sequence,
                                          campos[3], campos[4],
                                          campos[5], campos[6],
                                          campos[7])) # Instancia la clase 'Protein'.

    return(proteins)


def build_multifasta_for_muscle(estruct_dir, proteins, Protein_class):
    """
    Genera un archivo multifasta para cada query que contiene:
     - Los matches del BLAST
     - La seucencia query
    """
    # Crea el directorio donde se van a almacenar los multifasta.
    os.mkdir(estruct_dir.file_in_results_dir("unaligned_matches"))

    # Recorre la lista de queries con las que se ha instanciado la clase Protein
    for query in Protein_class.queries:
        # Abre el futuro fichero multifasta.
        with open(estruct_dir.file_in_results_dir("unaligned_matches/{}_matches.fa".format(query)), 'w') as out:
            for protein in proteins:
                # Identifica las instancias correspondientes a ese query.
                if protein.query_id == query:
                    out.write(">{}\n{}\n\n".format(protein.subject_id, protein.subject_seq))

            # Visita el fichero de queries para recuperar la secuencia query.
            with open(estruct_dir.file_in_data_dir("queries.fa"), 'r') as input_handle:
                for record in SeqIO.parse(input_handle, "fasta"):
                    if record.id == query:
                        sequence = str(record.seq)
                        header = record.id
                        out.write(">{}\n{}\n\n".format(header, sequence))
