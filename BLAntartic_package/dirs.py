#
# BLAntarctic v0.1 - 'dirs' module
#
# Este módulo contiene funciones y clases necesarias para la gestión de los
# directorios y ficheros de salida creados en la ejecución del programa.
#

import os
import re
from datetime import datetime


class EstructuraDirectorios:
    """
    Sus instancias contienen atributos y métodos para trabajar con el árbol
    de directorios que se crea durante la ejecución del programa.
    """

    def __init__(self):
        pass

    def add_results(self, relative_results_path):
        """
        Almacena el path de la carpeta de resutlados, relativo y absoluto.
        """
        self.rel_results_path = relative_results_path

    def add_data(self, relative_data_path):
        """
        Almacena el path de la carpeta de datos, relativo y absoluto.
        """
        self.rel_data_path = relative_data_path

    def file_in_results_dir(self, filename):
        """
        Devuelve el path relativo de un fichero de nombre 'filename' en la
        carpeta de resultados.
        """
        return "{}/{}".format(self.rel_results_path, filename)

    def file_in_data_dir(self, filename):
        """
        Devuelve el path relativo de un fichero de nombre 'filename' en la
        carpeta de datos.
        """
        return "{}/{}".format(self.rel_data_path, filename)



def calcula_indice(nombre_proyecto):
    """
    Comprueba si ya existen directorios con el nombre de proyecto y,
    si es así, calcula qué índice numérico hay que añadir para que,
    al crear los directorios necesarios, no se sobreescriba ninguno existente.
    """

    # Recupera una lista con los directorios y ficheros del cwd.
    files = os.listdir("..")
    directorios_involucrados = []

    # Busca aquellos elementos de la lista que contienen el nombre de proyecto.
    for file in files:
        if re.search(nombre_proyecto, file):
            m = re.search("calabaza", file)
            directorios_involucrados.append(file)

    # Extrae el índice numérico de esos directorios/ficheros que contienen el
    # nombre de proyecto.
    numeros = []
    for directorio in directorios_involucrados:
        try:
            numeros.append(int(directorio.split("_")[-1]))
        except:
            pass

    if len(directorios_involucrados) == 0:
        # Si no hay directorios con el nombre de proyecto, no hace falta índex.
        return(0)
    elif len(numeros) == 0:
        # Si el directorio que contiene el nombre de nombre_proyecto no tiene
        # índice numérico, el index será 1.
        return(1)
    else:
        # Si sí hay archivos con índice numérico, el índice del nuevo archivo
        # sera el mayor +1.
        maximo = max(numeros)
        return(maximo+1)


def crea_directorios(index, nombre_proyecto, *argv):
    """
    Crea los directorios con nombre basado en los 'argv*', con el índice numérico que se indique y el nombre de proyecto. Según la estructura:
                         'NombreProyecto_argv_index'
    Devuelve, además, una tupla con los paths relativos de los directorios que ha creado.
    """

    lista = []
    if index == 0:
        # Si el índice es 0, no se añade al nombre.
        for arg in argv:
            directorio = "../{}_{}".format(nombre_proyecto, arg)
            os.mkdir(directorio)
            lista.append(directorio)
    else:
        for arg in argv:
            directorio = "../{}_{}_{}".format(nombre_proyecto, arg, index)
            os.mkdir(directorio)
            lista.append(directorio)

    return(tuple(lista))
