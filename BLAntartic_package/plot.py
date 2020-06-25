# -*- coding: utf-8 -*-
#
# BLAntarctic v0.1 - 'plot' module
#
# Genera gráficos estáticos e interactivos.
#

import sys
import os
import re

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.image as mpimg
from matplotlib.widgets import Button
matplotlib.use('TkAgg')

from Bio import Seq
from Bio import SeqIO

from blastp import Protein
import prosite
import infowindow


def plot_domain(x1, x2, y1, text, ax, subject_id, hover_boxes):
    """
    Grafica un rectángulo coloreado, para representar un dominio conservado en
    una proteína.
    """
    height = 0.2
    width = x2-x1

    # Establece colores predefinidos para los dominios más frecuentes.
    color_code = {
            'PS00004':'limegreen',
            'PS00005':'limegreen',
            'PS00006':'limegreen',
            'PS00007':'limegreen',
            'PS00008':'gold',
    }

    try:
        color = color_code[text]
    except:
        color = "tab:orange"

    # Plotea el rectángulo
    rect1 = patches.Rectangle((x1, y1), width, height,
                              label = text + "&" + subject_id + "&" + str(x1) + "&" + str(x2),
                              color = color, linewidth=0, picker=5)

    if not hover_boxes == False:
        hover_boxes.append(ax.add_patch(rect1))
        # Lo añade a la lista que alimenta la función que lo dota de propiedades
        # interactivas.
    else:
        ax.add_patch(rect1)

    return(hover_boxes)


def plot_inter_domain(x1, x2, y1, ax, label, hover_boxes):
    """
    Grafica un rectángulo de color girs, para representar una región de una
    proteína subject que no coindice con ningún dominio conservado.
    """
    height = 0.2
    width = x2-x1
    rect1 = patches.Rectangle((x1, y1), width, height, linewidth=0, facecolor='grey', label=label)

    if not hover_boxes == False:
        hover_boxes.append(ax.add_patch(rect1))
    else:
        ax.add_patch(rect1)

    return(x1+width)


def plot_query(x1, x2, y1, ax):
    """
    Grafica un rectángulo de color azul, para representar una
    región de la proteína query que no coindice con ningún dominio conservado.
    """
    height = 0.2
    width = x2-x1
    rect1 = patches.Rectangle((x1, y1), width, height, linewidth=0, facecolor='tab:blue')
    ax.add_patch(rect1)

    return(x1+width)


def plot_spacer(x1, x2, y1, ax, color="grey"):
    """
    Grafica una línea gris, para representar espaciadores en el
    alineamiento múltiple.
    """
    width = x2-x1
    rect1 = patches.Rectangle((x1, y1+0.075 ), width, 0.05, linewidth=0, facecolor=color)
    ax.add_patch(rect1)

    return(x1+width)


def search_aligned_domains(dict_pattern, proteins, queries, threshold=8):
    """
    Busca dominios conservados en la secuencia ALINEADA de cada una de las
    proteínas (es decir, secuencias que incluyen gaps resultantes del
    alineamiento). De esta manera, se calcula la posición de dichos dominios
    dentro de la representación que se va a hacer.

    El parámetro 'threshold' permite excluir dominios muy pequeños, para que
    el gráfico no esté demasiado recargado.
    """

    # Para los subjects
    for protein in proteins:
        for pattern in list(dict_pattern.keys()):
            matches = re.finditer(pattern, str(protein.aligned_seq))
            for match in matches:
                match_length = match.end() - match.start()
                if match_length > threshold:
                    # Almacena resultado en la misma instancia de donde
                    # se ha obtenido la secuencia.
                    protein.add_aligned_domain(
                            dict_pattern[pattern],
                            match.start(),
                            match.end()
                            )

    # Para el query, ídem.
    for protein in queries:
        for pattern in list(dict_pattern.keys()):
            matches = re.finditer(pattern, str(protein.aligned_seq))
            for match in matches:
                match_length = match.end() - match.start()
                if match_length > threshold:
                    protein.add_aligned_domain(
                            dict_pattern[pattern],
                            match.start(),
                            match.end()
                            )


def analyze_alignment(proteins):
    """
    Determina qué regiones de la secuencia alineada de cada proteína son
    aminoácidos y qué regiones son gaps del alineamiento (guiones).
    """
    for protein in proteins:
        matches = re.finditer(r"[A-Z]+", str(protein.aligned_seq))
        for match in matches:
            protein.ilustrate_alignment(match.start(), match.end())


def static_plot(estruct_dir, proteins, queries, query):
    """
    Genera un plot estático, que se almacenará en la carpeta de resultados como
    PNG, y que se utilizará como miniatura para el menú interactivo.
    """

    # Tamaño figura y múltiples subplots
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_axes([0.2, 0, 0.7, 0.9])
    ax_names = fig.add_axes([0,0,0.2,1]) # Subplot a la izda., con los nombres

    y = 0
    DEF_QUERY = query
    protein_labels = [] # Almacena las etiquetas. Después se colocarán a la izquierda de cada proteína en orden.

    for protein in proteins:
        if protein.query_id == DEF_QUERY:
            protein_labels.append(( # Almacena etiqueta.
                    y,
                    protein.subject_accs,
                    protein.subject_spec))

            # Límite izquierdo y derecho de la proteína.
            x_min = protein.non_gapped[0][0]
            x_max = protein.non_gapped[-1][1]

            # Una línea a lo largo de toda la proteína. Después, muchas regiones
            # de la misma quedarán cubierta por 'plot_inter_domain'.
            plot_spacer(x_min, x_max, y, ax)

            # Rectángulos para las Non-Gapped-Regions (NGRs).
            for NGR in protein.non_gapped:
                plot_inter_domain(NGR[0], NGR[1], y, ax, protein.subject_id, False)

            # Por último, marca los dominios.
            for domain in protein.aligned_domains:
                plot_domain(domain[1], domain[2], y, domain[0], ax, protein.subject_id, False)

            y += 0.4 # Shift hacia arriba para la siguiente proteína.


    # Plotea el query correspondiente a esa representación.
    for query in queries:
        if query.id == DEF_QUERY:
            x_min = query.non_gapped[0][0]
            x_max = query.non_gapped[-1][1]
            plot_spacer(x_min, x_max, y, ax, color="tab:blue")
            for NGR in query.non_gapped:
                plot_query(NGR[0], NGR[1], y, ax)
            for domain in query.aligned_domains:
                plot_domain(domain[1], domain[2], y, domain[0], ax, protein.subject_id, False)
            ax.text(-10,y+0.085, "QUERY", ha='right', size='medium', va = 'center')

    # Coloca etiquetas a la izquierda de cada proteína, identificándola.
    for label in protein_labels:
        ax.text(-10,label[0]+0.07, label[1], ha='right', size='large', va = 'center')
        ax.text(-10,label[0]+0.13, label[2], ha='right', size='large', va = 'center')

    ax.axis('off')  # Quita ejes
    ax_names.axis("off")  # Quita ejes
    ax.autoscale()
    ax.set_title(str(DEF_QUERY))

    # Calcula el tamaño de la figura según el número de proteínas que contiene.
    fig_width = len(protein_labels)
    fig.set_size_inches(fig_width*1.33, fig_width, forward=True)

    # Guarda la figura.
    path = estruct_dir.file_in_results_dir("plots/{}.png".format(DEF_QUERY))
    plt.savefig(path)
    plt.close()


def generate_static_graphs(estruct_dir, proteins, queries, Protein_class):
    """
    Llama a la función anterior para cada query, y almacena los resultados en
    la carpeta correspondiente.
    """
    os.mkdir(estruct_dir.file_in_results_dir("plots"))
    for query in Protein_class.queries:
        static_plot(estruct_dir, proteins, queries, query)


def interactive_plot(selection, proteins, queries, dict_names, estruct_dir):
    """
    Genera un plota análogo al de la función 'static_plot', pero con un
    comportamiento interactivo.

    Se han omitido los comentarios en las partes que son idénticas a la función
    'static_plot'.
    """
    fig = plt.figure(1)
    ax = fig.add_axes([0.2, 0, 0.7, 0.9])
    ax_names = fig.add_axes([0,0,0.2,1])

    y = 0 # Init
    hover_boxes = []
    protein_labels = []
    DEF_QUERY = selection

    for protein in proteins:
        if protein.query_id == DEF_QUERY:
            protein_labels.append((
                    y,
                    protein.subject_accs,
                    protein.subject_spec))

            x_min = protein.non_gapped[0][0]
            x_max = protein.non_gapped[-1][1]
            plot_spacer(x_min, x_max, y, ax)

            for NGR in protein.non_gapped:
                plot_inter_domain(NGR[0], NGR[1], y, ax, protein.subject_id, False)

            for domain in protein.aligned_domains:
                plot_domain(domain[1], domain[2], y, domain[0], ax, protein.subject_id, hover_boxes)
                # El parámetro 'hover_boxes' indica que se debe añadir el
                # rectángulo a la lista sobre la cual trabaja la función
                # 'hover'.

            y += 0.4

    for query in queries:
        if query.id == DEF_QUERY:
            x_min = query.non_gapped[0][0]
            x_max = query.non_gapped[-1][1]
            plot_spacer(x_min, x_max, y, ax, color="tab:blue")
            for NGR in query.non_gapped:
                plot_query(NGR[0], NGR[1], y, ax)

            for domain in query.aligned_domains:
                plot_domain(domain[1], domain[2], y, domain[0], ax, protein.subject_id, hover_boxes)
            ax.text(-10,y+0.085, "QUERY", ha='right', size='medium', va = 'center')

    # Coloca etiquetas a la izquierda de cada proteína, identificándola.
    for label in protein_labels:
        ax.text(-10,label[0]+0.05, label[1], ha='right', size='x-small', va = 'center')
        ax.text(-10,label[0]+0.15, label[2], ha='right', size='x-small', va = 'center')


    # Define cómo será la anotación que se muestra al hacer hover con el ratón
    # sobre un dominio.
    annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    def update_annot(artist):
        """
        Actualiza el contenido de la anotación, según sobre qué dominio esté
        posado el ratón.
        """
        center_x = artist.get_x() + artist.get_width() / 2
        center_y = artist.get_y() + artist.get_height() / 2
        annot.xy = (center_x, center_y)

        retrieved_label = str(artist.get_label())
        split_label = retrieved_label.split("&")[0]

        if "@" in split_label:
            text = split_label
        else:
            text = dict_names[split_label] + " · " + split_label

        annot.set_text(text)
        annot.get_bbox_patch().set_alpha(0.4)


    def hover(event):
        """
        Cuando el ratón está sobre un dominio, muestra la anotación, y llama
        a la función 'update_annot' para actualizar su contenido.
        """
        vis = annot.get_visible()
        if event.inaxes == ax:
            an_artist_is_hovered = False
            for artist in hover_boxes: # Considera los dominios de la lista.
                contains, _ = artist.contains(event)
                if contains:
                    an_artist_is_hovered = True
                    update_annot(artist)
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
            if not an_artist_is_hovered:
                # Si el ratón no está sobre un dominio, oculta la anotación.
                annot.set_visible(False)
                fig.canvas.draw_idle()

    def on_pick(event):
        """
        Abre el pop-up de información extendida si se hace clic sobre un
        dominio.
        """
        print("Ha hecho clic en un dominio. Se está abriendo la ventana de información expandida. (Esto puede tardar unos segundos.)")

        selection = event.artist.get_label()
        fields = selection.split("&")

        infowindow.call_expanded_info(fields[0], fields[1], (fields[2], fields[3]), estruct_dir) # Llama a la función que abre el pop-up.

    def volver(event):
        """
        Acción del botón 'volver', para regresar al menú de miniaturas.
        """
        plt.close() # Cierra el plot actual.
        print("Volviendo al menú...")
        # Abre el menú de nuevo.
        thumbnails_menu(estruct_dir, proteins, queries, dict_names)

    # Captura la actividad del ratón sobore los dominios
    fig.canvas.mpl_connect("motion_notify_event", hover) # Posarse encima
    fig.canvas.callbacks.connect('pick_event', on_pick)  # Clic

    # Botón volver atrás
    axcut = plt.axes([0.9, 0.0, 0.1, 0.075])
    bcut = Button(axcut, 'Volver', color='tab:grey', hovercolor='silver')
    bcut.on_clicked(volver)

    ax.axis('off')
    ax_names.axis("off")
    ax.autoscale()
    ax.set_title(str(DEF_QUERY))

    plt.show()


def calculate_geometry(n):
    """
    Calcula la disposición de subplots necesaria para acomodar tantas miniaturas
    como deba contener el menú.
    """
    if n > 4:
        cols = 4
        rows = int(n/cols)
        if n%cols > 0:
            rows += 1

        enteras = int(n/4)*4
        odds = n - enteras

        plotting_range = [i+1 for i in range(enteras)]

        if odds == 1:
            plotting_range.append(plotting_range[-1]+2)
        elif odds == 2:
            plotting_range.append(plotting_range[-1]+2)
            plotting_range.append(plotting_range[-1]+1)
        elif odds == 3:
            for i in range(3):
                plotting_range.append(plotting_range[-1]+1)

        return(rows, plotting_range)

    else: # Si hay menos de cuatro elementos, una única fila.
        plotting_range = [i+1 for i in range(n)]
        return(1, plotting_range)


def thumbnails_menu(estruct_dir, proteins, queries, dict_names):
    """
    Genera un menú con miniatura de los gráficos disponibles, de forma que el
    usuario pueda hacer clic sobre uno de ellos y abrirlo en su versión
    interactiva.
    """
    # Recupera los paths de las miniaturas (PNGs).
    path = estruct_dir.file_in_results_dir("plots")
    thumbs = []
    for filename in os.listdir(path):
        thumbs.append(path+"/"+filename)

    # Establece la disposición de subplots.
    n = len(thumbs)
    rows, plotting_range = calculate_geometry(n)
    plot_height = rows * 2
    fig = plt.figure(figsize=(10,plot_height))

    # Imprime en cada subplot la imagen correspondiente.
    for i, path in zip(plotting_range, thumbs):
        # Crea subplot.
        plt.subplot(rows, 4, i)
        ax = plt.gca()
        ax.axis("off")

        # Inserta imagen.
        name = path.split("/")[-1]
        name = name.split(".")[0]
        image = mpimg.imread(path)
        ax.imshow(image, label = name, picker=5)

        # Crea título para la miniatura.
        xlim = ax.get_xlim()
        ax.text(xlim[1]/2, -1, name)

    plt.margins(4)
    plt.suptitle("Haz clic en un alineamiento para abrirlo.")


    def on_pick(event):
        """
        Abre el gráfico interactivo correspondiente cuando se hace clic en la
        miniatura.
        """
        plt.close() # Cierra menú.
        selection = event.artist.get_label()
        # Abre el plot interactivo.
        interactive_plot(selection, proteins, queries, dict_names, estruct_dir)

    # Captura clics del ratón.
    fig.canvas.callbacks.connect('pick_event', on_pick)

    plt.show()
