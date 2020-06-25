# -*- coding: utf-8 -*-
#
# BLAntarctic v0.1 - 'prosite' module
#
# Búsqueda de dominios conservados.
#

import os
import re

from Bio.ExPASy import Prosite
from Bio.ExPASy import Prodoc

def convert_RE(expression):
    """
    Convierte las expresiones regulares del fichero 'prosite.dat' al formato
    reconocido por el paquete 're' de Python.
    """
    dict = {
        '-' : '',
        '{' : '[^',
        '}' : ']',
        '(' : '{',
        ')' : '}',
        'x' : '.',
        '>' : '$',
        '<' : '^'
        }

    for original, translation in dict.items():
        expression = expression.replace(original, translation)

    return str(expression)

def create_prosite_dict(exclude = False):
    """
    Crea dos diccionarios: uno que relaciona los patrones RE de cada dominio
    con el número de accessión del mismo ('dict_pattern'), y otro que relaciona el número de accesión con el nombre común ('dict_names').

    Además, si el parámetro 'exclude' es 'True', se excluyen del diciconario
    dominios con alta probabilidad de ocurrencia (sistios de fosforilación,
    etc.), que contienen la etiqueta '/SKIP-FLAG=TRUE>' (ver manual de ScanProsite para más información.)
    """
    dict_pattern = {}
    dict_names = {}
    with open('./prosite.dat','r') as handle:
        records = Prosite.parse(handle)
        for record in records:
            if record.pattern != "": # Hay algunos records que no tienen patrón.
                if exclude == False:
                    converted_RE = convert_RE(record.pattern)
                    # Recupera el nombre y accessión del dominio.
                    dict_pattern[converted_RE] = str(record.accession)
                    dict_names[str(record.accession)] = str(record.name)
                else:
                    if not record.cc_skip_flag == "TRUE":
                        converted_RE = convert_RE(record.pattern)
                        # Recupera el nombre y accessión del dominio.
                        dict_pattern[converted_RE] = str(record.accession)
                        dict_names[str(record.accession)] = str(record.name)

    return dict_pattern, dict_names


def search_domains(estruct_dir, proteins, Protein_class, dict_pattern):
    """
    Recorre todas las instancias de 'Protein', y busca en su secuencia de
    aminoácidos coincidencias con alguna de las expresiones regulares de
    el diccionario de Prosite.

    Almacena los matches en la lista 'self.domains' de la instancia
    correspondiente, en forma de tuplas (accession, start, end).
    """
    for query in Protein_class.queries:
        for protein in proteins:
            if protein.query_id == query:
                for pattern in list(dict_pattern.keys()):
                    ungapped_seq = protein.subject_seq.replace("-", "")
                    matches = re.finditer(pattern, ungapped_seq)
                    for match in matches:
                        protein.add_domain(
                                dict_pattern[pattern],
                                match.start(),
                                match.end()
                                )


def output_domains(estruct_dir, proteins, Protein_class, dict_names):
    """
    Genera un archivo de texto separado con tabulaciones con los dominios
    encontrados para cada query.
    """
    os.mkdir(estruct_dir.file_in_results_dir("protein_domains"))

    for query in Protein_class.queries:
        with open(estruct_dir.file_in_results_dir("protein_domains/{}_domains.tsv".format(query)), 'w') as out:

            out_txt = ""
            out_txt += "Query\tProtein\tProsite_accs\tDomain_name\tSpan\n\n"
            for protein in proteins:
                if query == protein.query_id:
                    for domain in protein.domains:
                        accession = domain[0]
                        name = dict_names[accession]
                        out_txt += "{}\t{}\t{}\t{}\t{}-{}\n".format(
                                query,
                                protein.subject_id_spaced,
                                accession,
                                name,
                                domain[1],
                                domain[2]
                        )
                    out_txt += "\n"

            out.write(out_txt)
