# -*- coding: utf-8 -*-
# BLAntarctic v0.1 - 'infowindow' module
#
# Funciones para generar el pop-up informativo cuando se hace clic sobre un
# dominio en el plot interactivo.
#

import tkinter as tk
from tkinter import ttk
import os
import textwrap

from Bio.ExPASy import Prosite,Prodoc
from Bio import Seq
from Bio import SeqIO


def retrieve_form_prosite(accession):
    """
    Recupera de prosite.dat y prosite.doc la información relativa al dominio
    conservado que se va a mostrar en la ventana.
    """
    with open('./prosite.dat','r') as handle:
        dat_records = Prosite.parse(handle)
        for record in dat_records:
            if record.accession == accession:
                name = record.name
                description = record.description
                pattern = record.pattern
                doc_accession = record.pdoc

    with open('./prosite.doc','r',  encoding='latin1') as handle:
        doc_records = Prodoc.parse(handle)
        for record in doc_records:
            if record.accession == doc_accession:
                text = record.text

    return(name, description, pattern, text)


def retrieve_from_GB(estruct_dir, protein):
    """
    Recupera del fichero GeneBank correspondiente la información relativa a la
    proteína donde se encuentra el dominio y el genoma del que procede, para mostrarlo en la ventana.
    """

    accession, species = protein.split("@")

    path = estruct_dir.file_in_data_dir("raw_GBs")

    for filename in os.listdir(path):
        if filename == species+".gbff":
            with open("{}/{}".format(path, filename), "r") as input_handle:
                for record in SeqIO.parse(input_handle, "genbank"):
                    for feature in record.features:
                        if (feature.type == 'CDS'
                            and feature.qualifiers['locus_tag'][0] == accession):
                            location = feature.location

                            try:
                                EC_number = feature.qualifiers['EC_number'][0]
                            except:
                                EC_number = "N/A"

                            try:
                                product = feature.qualifiers['product'][0]
                            except:
                                product = "N/A"

                            try:
                                raw_translation = feature.qualifiers['translation'][0]
                                translation = textwrap.fill(raw_translation, 58, break_long_words=True)
                            except:
                                translation = "N/A"

    return(species, accession, location, EC_number, product, translation)



def call_expanded_info(prosite_accession, protein, span, estruct_dir):
    """
    Abre la ventana. Esta función es la que se llama al hacer click sobre un
    dominio en la vista interactiva.
    """
    main_window = tk.Tk()
    app = Application(main_window, str(prosite_accession), span, estruct_dir, protein)
    app.mainloop()


class DomainFrame(ttk.Frame):
    """
    Tikinter. Pestaña con la información del dominio.
    """

    def __init__(self, notebook, prosite_accession, span):
        super().__init__(notebook)

        # Llama a la función retrieve_from_prosite para recuperar la información
        # que apareecerá en la pestaña.
        (name, description,
         pattern, text) = retrieve_form_prosite(prosite_accession)

        position = "{} - {}".format(span[0],span[1])

        # ----- doc_accession
        self.accs_title = ttk.Label(self, text = "Accession:",
                                     font='Helvetica 14 bold')
        self.accs_title.pack(anchor='nw')

        self.accs_number = ttk.Label(self, text = prosite_accession)
        self.accs_number.pack(anchor='nw')

        self.sep = ttk.Separator(self, orient="horizontal")
        self.sep.pack(anchor="nw", fill=tk.X)


        # ----- Name
        self.name_title = ttk.Label(self, text = "Name:",
                                     font='Helvetica 14 bold')
        self.name_title.pack(anchor='nw')

        self.name_number = ttk.Label(self, text = name)
        self.name_number.pack(anchor='nw')

        self.sep = ttk.Separator(self, orient="horizontal")
        self.sep.pack(anchor="nw", fill=tk.X)


        # ----- Description
        self.description_title = ttk.Label(self, text = "Description:",
                                     font='Helvetica 14 bold')
        self.description_title.pack(anchor='nw')

        self.description_number = ttk.Label(self, text = description)
        self.description_number.pack(anchor='nw')

        self.sep = ttk.Separator(self, orient="horizontal")
        self.sep.pack(anchor="nw", fill=tk.X)

        # ----- Pattern
        self.pattern_title = ttk.Label(self, text = "Pattern:",
                                     font='Helvetica 14 bold')
        self.pattern_title.pack(anchor='nw')

        self.pattern_number = ttk.Label(self, text = pattern)
        self.pattern_number.pack(anchor='nw')

        self.sep = ttk.Separator(self, orient="horizontal")
        self.sep.pack(anchor="nw", fill=tk.X)

        # ----- Pattern
        self.pattern_title = ttk.Label(self, text = "Position in alignment:",
                                     font='Helvetica 14 bold')
        self.pattern_title.pack(anchor='nw')

        self.pattern_number = ttk.Label(self, text = position)
        self.pattern_number.pack(anchor='nw')

        self.sep = ttk.Separator(self, orient="horizontal")
        self.sep.pack(anchor="nw", fill=tk.X)

        # ----- Text
        self.pattern_number = ttk.Label(self, text = text)
        self.pattern_number.pack(anchor='nw')


class ProteinFrame(ttk.Frame):
    """
    Tkinter. Pestaña con la información de la proteína en su contexto.
    """

    def __init__(self, notebook, estruct_dir, protein):
        super().__init__(notebook)

        # Llama a la función retrieve_from_GB para recuperar la información
        # que apareecerá en la pestaña.
        (species, accession,
         location, EC_number,
         product, translation) = retrieve_from_GB(estruct_dir, protein)

        # ----- Species
        self.species_title = ttk.Label(self, text = "Species:",
                                     font='Helvetica 14 bold')
        self.species_title.pack(anchor='nw')

        self.species_number = ttk.Label(self, text = species)
        self.species_number.pack(anchor='nw')

        self.species = ttk.Separator(self, orient="horizontal")
        self.species.pack(anchor="nw", fill=tk.X)

        # ----- Locus-tag
        self.locus_title = ttk.Label(self, text = "Locus tag:",
                                     font='Helvetica 14 bold')
        self.locus_title.pack(anchor='nw')

        self.locus_number = ttk.Label(self, text = accession)
        self.locus_number.pack(anchor='nw')

        self.locus = ttk.Separator(self, orient="horizontal")
        self.locus.pack(anchor="nw", fill=tk.X)

        # ----- Location
        self.location_title = ttk.Label(self, text = "Location:",
                                     font='Helvetica 14 bold')
        self.location_title.pack(anchor='nw')

        self.location_number = ttk.Label(self, text = location)
        self.location_number.pack(anchor='nw')

        self.location = ttk.Separator(self, orient="horizontal")
        self.location.pack(anchor="nw", fill=tk.X)

        # ----- EC_number
        self.ec_title = ttk.Label(self, text = "EC number:",
                                     font='Helvetica 14 bold')
        self.ec_title.pack(anchor='nw')

        self.ec_number = ttk.Label(self, text = EC_number)
        self.ec_number.pack(anchor='nw')

        self.ec = ttk.Separator(self, orient="horizontal")
        self.ec.pack(anchor="nw", fill=tk.X)

        # ----- Product
        self.product_title = ttk.Label(self, text = "Product:",
                                     font='Helvetica 14 bold')
        self.product_title.pack(anchor='nw')

        self.product_number = ttk.Label(self, text = product)
        self.product_number.pack(anchor='nw')

        self.product = ttk.Separator(self, orient="horizontal")
        self.product.pack(anchor="nw", fill=tk.X)


        # ----- Translation
        self.translation_title = ttk.Label(self, text = "Translation:",
                                     font='Helvetica 14 bold')
        self.translation_title.pack(anchor='nw')

        self.translation_number = ttk.Label(self, text = translation)
        self.translation_number.pack(anchor='nw')

        self.translation = ttk.Separator(self, orient="horizontal")
        self.translation.pack(anchor="nw", fill=tk.X)

class Application(ttk.Frame):
    """
    Tkinter. Genera la ventana (notebook) que contiene las dos pestañas
    definidas en las clases anteriores.
    """

    def __init__(self, main_window, prosite_accession, span, estruct_dir, protein):
        super().__init__(main_window)
        main_window.title("Expanded info - BLAntarctic v0.1")

        self.notebook = ttk.Notebook(self)

        self.domain_frame = DomainFrame(self.notebook, prosite_accession, span)
        self.notebook.add(
            self.domain_frame, text="Domain", padding=10)

        self.protein_frame = ProteinFrame(self.notebook, estruct_dir, protein)
        self.notebook.add(
            self.protein_frame, text="Protein", padding=10)

        self.notebook.pack(padx=10, pady=10)
        self.pack()
