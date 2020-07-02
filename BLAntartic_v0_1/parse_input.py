# -*- coding: utf-8 -*-
#
# BLAntarctic v0.1 - 'parse_input' module
#
# Procesa los ficheros input (queries, subjects).-
#

import sys
import os
import shutil

import re
from Bio import Seq
from Bio import SeqIO


def make_genomes_multifasta(folder_path, estruct_dir):
    """
    Parsea todos los GBs contenido en la carpeta 'folder_path', extrae las
    secuencias de proteína de sus CDSs y econstruye un multifasta conjunto.

    En el header de cada secuencia del MULTIFASTA, incluye el nombre del gen,
    así como la especie del GB de donde se ha estraído (separado por '@').
    """
    multifasta_path = estruct_dir.file_in_data_dir('genomes_multifasta.fa')
    multifasta = open(multifasta_path, 'a')

    destination_path = estruct_dir.file_in_data_dir("raw_GBs/")
    os.mkdir(destination_path)

    genomes_dict = {}

    for filename in os.listdir(folder_path):
        try:
            with open("{}/{}".format(folder_path, filename), "r") as input_handle:

                for record in SeqIO.parse(input_handle, "genbank"):
                    # Extrae el nombre de la especie (bionmial) y sustituye
                    # espacios por guiones bajos.
                    nombre_largo = record.description.split()
                    nombre_corto = nombre_largo[0] + "_" + nombre_largo[1]

                    # Construye el MULTIFASTA.
                    for feature in record.features:
                        if feature.type == 'CDS':
                            try:
                                multifasta.write(">{}@{}\n{}\n\n".format(feature.qualifiers['locus_tag'][0],nombre_corto, feature.qualifiers['translation'][0]))
                            except:
                                pass
                shutil.copy2("{}/{}".format(folder_path, filename),
                              destination_path+"/"+nombre_corto+".gbff")
        except:
            # Si el parser de BioPy no puede abrir el archivo...
            pass
    multifasta.close()


def store_GBs_copies(folder_path, estruct_dir):
    """
    Crea una copia de los GeneBanks originales en la carpeta 'data/raw_GBs',
    para que estén accesibles si el programa se vuelve a ejecutar a partir de
    un fichero .bapj (BLAntartic Project). Estos ficheros sirven de input para
    la pestaña "Protein" de la ventana de información expandida que se abre al
    hacer clic sobre un dominio proteico en la vista interactiva.
    """

    destination_path = estruct_dir.file_in_data_dir("raw_GBs/")
    os.mkdir(destination_path)

    for filename in os.listdir(folder_path):
        origin = folder_path + "/" + filename
        shutil.copy2(origin, destination_path)

def preprocess_queries(query_multifasta, estruct_dir):
    """
    Recorta el header de los queries a un máximo de 6 caracteres y añade un
    índice numérico ('Q1, Q2, Q3...') para diferenciarlos en el caso de que
    los 6 primeros caracteres coincidan en varios queries.

    Crea un farchivo multifasta con estas modificaciones en la carpeta de
    'data'.
    """

    out_path = estruct_dir.file_in_data_dir("queries.fa")
    out = open(out_path, 'w')

    queries = []

    with open(query_multifasta, "r") as input_handle:
        record_index = 1
        for record in SeqIO.parse(input_handle, "fasta"):
            sequence = str(record.seq)
            header = record.id.replace("_","")

            regex = re.compile('[^a-zA-Z]')
            clean_header = regex.sub('', header)

            if len(clean_header) > 6:
                clean_header = clean_header[0:6]
            header1 = ">Q"+str(record_index) +"_"+ clean_header
            out.write(header1+"\n")
            out.write(sequence+"\n\n")


            record_index += 1
    out.close()

    return(queries)


def check_dbs():
    """
    Comprueba que los ficheros de la base de datos PROSITE están contenidos en
    el directorio de trabajo.
    """

    if not "prosite.dat" in os.listdir() or not "prosite.doc" in os.listdir():
        print("ERROR. En la carpeta que contiene los scripts (directorio de trabajo) se deben incluir los ficheros 'prosite.dat' y 'prosite.doc'")
        sys.exit(2)
