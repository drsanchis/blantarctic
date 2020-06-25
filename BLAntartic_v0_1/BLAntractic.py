#!/usr/bin/env python3
#
# BLAntarctic v0.1
# Diego Ruiz Sanchis, 2020
# Programación para la bioionformática
# 

import getopt
import sys
import os

import blastp
import dirs
import muscle
import parse_input
import prosite
import plot
import saveproject


script_name = sys.argv[0]
argv0 = sys.argv[1:]

def use(script_name):
    """
    Muestra una ayuda breve del script.
    """

    print("""USO:
          {} -n name -s subjects_directory -q queries.multifasta

            · 'name' es el nombre de proyecto.

            · 'subjects_directory' es un directorio que contiene los ficheros
              '.gbff ' de los genomas que se quieren usar como subjects.

            · 'queries.multifasta' es un fichero (multi)fasta que contiene
              la(s) secuencia(s) de la(s) proteína(s) que se quiere(n)
              utilizar como query/queries.

            opcional: [-h] [--eval blast_evalue] [--ident blast_identity]
                    [--cov blast_coverage] [--exclude] [-o project.bapj]""".format(script_name))


def help():
    """
    Imprime en pantalla la ayuda.
    """
    with open("./help.txt", 'r') as filehandle:
        print(filehandle.read())
    sys.exit(0)


# Valores por defecto si no se cambian a través de las opciones.
eval_threshold = 0.00001
cov_threshold = 50
ident_threshold = 30
exclude = False
opening = False

# Extrae las opciones y argumentos posicionales.
try:
    opts, argumentos = getopt.getopt(argv0, 'n:hs:q:o:', ['eval=', 'ident=', 'cov=', 'exclude'])

except getopt.GetoptError:
    print("ERROR: Error en las opciones o argumentos introducidos.")
    use(script_name)
    sys.exit(2)

for o, a in opts:
    if o in ("-h", "--help", "-help"):
        help()

    elif o == "-s":
        subjects_directory = str(a)

    elif o == "-q":
        queries_file = str(a)

    elif o == "-n":
        project_name = str(a)

    elif o == "-o":
        opening = True
        open_path = str(a)

    elif o == "--evalue":
        eval_threshold = float(a)

    elif o == "--ident":
        ident_threshold = float(a)

    elif o == "--cov":
        cov_threshold = float(a)

    elif o in "--exclude":
        exclude = True

if opening == False:
    try:
        subjects_directory
        queries_file
        project_name
    except:
        print("ERROR: Los argumentos '-s', '-q' y '-n' son necesarios.")
        use(script_name)
        sys.exit(2)


# Comprueba que la base de datos de prosite está contenida en el wd.
parse_input.check_dbs()

if opening == False:

    # Crea los directorios de resultados y datos.
    print("Creando directorios...")
    indice = dirs.calcula_indice(project_name)
    path_results, path_data = dirs.crea_directorios(indice, project_name,  "results", "data")

    # Almacena los paths de los directorios creados en una instancia de la clase
    # 'EstructuraDirectorios'.
    estruct_dir = dirs.EstructuraDirectorios()
    estruct_dir.add_results(path_results)
    estruct_dir.add_data(path_data)

    # Parsea los fichero de input y extrae secuencias, etc.
    print("Procesando inputs...")
    parse_input.make_genomes_multifasta(subjects_directory, estruct_dir)
    queries = parse_input.preprocess_queries(queries_file, estruct_dir)
    # parse_input.store_GBs_copies(subjects_directory, estruct_dir)

    # Hace el BLAST y almacena los resultados en un fichero XML.
    print("Haciendo BLAST...")
    blastp.blast_it(estruct_dir, eval_threshold, cov_threshold)

    # Filtra los resultados del blast por identidad.
    print("Filtrando BLAST...")
    blastp.filter_ident(estruct_dir, ident_threshold)

    # Almacena cada resultado del blast en una instancia de la clase Protein.
    print("Instanciando resultados de BLAST...")
    Protein_class = blastp.Protein
    blast_results = blastp.parse_blast_results(estruct_dir, Protein_class)

    # Crea un fichero multifasta con los resultados del blast para cada query.
    print("Construyendo multifasta para MUSCLE...")
    blastp.build_multifasta_for_muscle(estruct_dir, blast_results, Protein_class)

    # Hace un alineamiento múltiple sobre cada uno de los multifasta del paso
    # anterior.
    print("Haciendo alineamiento múltiple...")
    muscle.multi_align(estruct_dir)
    aligned_queries = muscle.parse_alignment(estruct_dir, blast_results)


    # Construye árbol filogenético en formato newick.
    print("Construyendo árbol filogenético...")
    muscle.build_tree(estruct_dir)

    # Representa cada árbol filogenético en un fichero .pdf.
    muscle.plot_tree(estruct_dir)

    # Parsea la base de datos de prosite (prosite.dat) y genera dos diccionarios que
    # relacionan los patrones con los números de accesión y nombre de los dominios.
    print("Parseando base de datos PROSITE...")
    dict_pattern, dict_names = prosite.create_prosite_dict(exclude=exclude)

    # Busca dominios conservados en los resultados del blast.
    print("Buscando dominios conservados...")
    prosite.search_domains(estruct_dir, blast_results, Protein_class, dict_pattern)

    # Genera un fichero output para cada query con un resumen de los dominios
    # encontrados
    prosite.output_domains(estruct_dir, blast_results, Protein_class, dict_names)


    print("Ubicando dominios conservados en el alineamiento...")
    plot.search_aligned_domains(dict_pattern, blast_results, aligned_queries, threshold=4)

    print("Construyendo plot...")
    plot.analyze_alignment(blast_results)
    plot.analyze_alignment(aligned_queries)

    print("Guardando proyecto...")
    saveproject.save(estruct_dir, blast_results, aligned_queries, Protein_class, dict_names)

    print("Generando plots...")
    plot.generate_static_graphs(estruct_dir, blast_results, aligned_queries, Protein_class)

if opening == True:
    Protein_class = blastp.Protein
    objects = saveproject.open_project(open_path)
    estruct_dir = objects[0]
    blast_results = objects[1]
    aligned_queries = objects[2]
    Protein_class.import_queries_list(objects[3])
    dict_names = objects[4]


print("Abriendo menú interactivo...")
plot.thumbnails_menu(estruct_dir, blast_results, aligned_queries, dict_names)
