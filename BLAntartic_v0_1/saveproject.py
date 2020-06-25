# -*- coding: utf-8 -*-
#
# BLAntarctic v0.1 - 'saveproject' module
#
# Gestión de ficheros *.bapj.
#

import pickle

def save(estruct_dir, proteins, queries, Protein_class, dict_names):
    """
    Vuelca en un fichero el contenido de los objectos necesarios para reanudar
    la ejecución del script más tarde.
    """
    path = estruct_dir.file_in_results_dir("project.bapj")
    queries_list = Protein_class.queries
    all_objects = (estruct_dir, proteins, queries, queries_list, dict_names)
    with open(path, 'wb') as output:
        pickle.dump(all_objects, output, pickle.HIGHEST_PROTOCOL)

def open_project(path):
    """
    Recupera de un fichero los objetos necesarios para reanudar la ejecución.
    """
    with open(path, 'rb') as input:
        all_objects = pickle.load(input)

    return(all_objects)
