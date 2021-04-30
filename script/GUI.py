#!/usr/bin/env python3
from tkinter import *
from tkinter import ttk

import os
import time

import pandas as pd

import fetch as fetch

class GUI:
    def __init__(self):
        # attribute
        self.organism_df = fetch.load_df_from_pickle()
        self.is_in_critical_section = False
        # Creating window
        self.window = Tk()
        self.window.geometry("1200x800")
        self.window.title("GenBank Application")

        # Frame 1 : Liste des fichiers
        (self.treeview, self.scrollbar) = self.create_tree()
        self.treeview.bind("<<TreeviewSelect>>", self.on_tree_select)

        # Frame 2 : Informations
        (self.Info, self.Info_lab, self.labelText, self.depositLabel) = self.create_info()

        # Frame 3 : Logs
        (self.Logs, self.Logs_lab, self.log_text, self.log_scroll) = self.create_log()

        # Regions menu
        (self.selected_region, self.menu_menu, self.run_search, self.OptionList_region) = self.create_region_menu()

        # Reset button
        self.reset_button = self.create_reset_button()

        # Start mainloop
        self.window.mainloop()

    # TREE CREATION METHODS
    def tree_exist(self, tree, name) :
        if tree.exists(name) :
            name += '1'
            self.tree_exist(tree, name)
        return name

    def create_node(self, tree, node_path, node) :
        name_nodes = os.listdir(node_path)
        if name_nodes == [] : pass
        for name in name_nodes :
            name_path = node_path + name
            if os.path.isdir(name_path) :
                name = self.tree_exist(tree, name)
                name_path += '/'
                try :
                    tree.insert(node, '1000000',iid=name, text=name)
                    self.create_node(tree, name_path, name)
                except:
                    print(name_path)
                    pass


    def create_tree(self):
        list = Frame(self.window)
        list.place(x=0, y=0, anchor="nw", width=400, height=800)

        # Creating treeview window
        scrollbar = Scrollbar(list)
        scrollbar.pack( side = RIGHT, fill = Y )

        treeview = ttk.Treeview(list)
        treeview.heading('#0', text='Arborescence des fichiers')

        treeview.configure(yscrollcommand=scrollbar.set)
        treeview.pack(fill="both", expand=True)
        scrollbar.configure(command=treeview.yview)

        root_path = "../Results/"
        treeview.insert('', '0', text='Results', iid='Results')
        self.create_node(treeview, root_path, 'Results')

        return treeview, scrollbar

    # INFO CREATION METHODS
    def create_info(self):
        Info = Frame(self.window)
        Info.place(x=400, y=0, anchor="nw", width=800, height=600)

        Info_lab = LabelFrame(Info, text="Informations", padx=20, pady=20)
        Info_lab.pack(fill="both", expand="yes")

        Label(Info_lab, text="Organisme choisi : ").grid(row = 0, column = 0, sticky = W, ipadx = 100, ipady = 30)

        labelText = StringVar()
        labelText.set('Aucun')

        depositLabel = Label(Info_lab, textvariable=labelText)
        depositLabel.grid(row = 0, column = 1, sticky = W)

        Label(Info_lab, justify = LEFT, text="Région fonctionnelle choisie : ").grid(row = 4, column = 0, sticky = W, ipadx = 100, ipady = 30)
        return Info, Info_lab, labelText, depositLabel

    # LOG CREATION METHODS
    def create_log(self):
        Logs = Frame(self.window, background="#b22222")
        Logs.place(x=400, y=400, anchor="nw", width=800, height=400)

        Logs_lab = LabelFrame(Logs, text="Logs")
        Logs_lab.pack(fill="both", expand="yes")


        log_text = Text(Logs_lab, height='400', width='800')
        log_scroll = Scrollbar(Logs_lab, command = log_text.yview)
        log_text.configure(yscrollcommand=log_scroll.set)

        log_text.pack(side=LEFT)
        log_scroll.pack(side=RIGHT, fill=Y)
        return Logs, Logs_lab, log_text, log_scroll

    # REGION MENU CREATION METHODS
    def create_region_menu(self):
        OptionList = [
        "Aucun",
        "CDS",
        "centromere",
        "intron",
        "mobile_element",
        "ncRNA",
        "rRNA",
        "telomere",
        "tRNA",
        "3'UTR",
        "5'UTR"
        ]

        variable = StringVar(self.Info_lab)
        variable.set(OptionList[0])

        menu = OptionMenu(self.Info_lab, variable, *OptionList)
        menu["borderwidth"]=1
        menu.grid(row = 4, column = 1, sticky = W)

        variable.trace("w", self.callback)

        run_search = Button(self.Info_lab, text ="Search", command = self.search_button_callback, relief = RIDGE, borderwidth=1)
        run_search.grid(row = 5, sticky = 'se', column = 2, ipadx = 20, pady = 30, padx = 30)
        return variable, menu, run_search, OptionList

    # RESET BUTTON CREATION
    def create_reset_button(self):
        reset_button = Button(self.Info_lab, text ="Reset Tree", command = self.reset_button_callback, relief = RIDGE, borderwidth=1)
        reset_button.grid(row = 7, sticky = 'se', column = 2, ipadx = 20, pady = 30, padx = 30)
        return reset_button

    # FEATURES METHODS
    def print_on_window(self, t): #affiche t dans les logs
        time_string = time.strftime('%H:%M:%S')
        self.log_text.insert(INSERT, time_string + ' : ' + t + "\n")
        self.log_text.yview(END)

    def callback(self, *args): # fonction pour executer du code pour le menu, a changer
        self.print_on_window("The selected item is " + self.selected_region.get())

    def get_path(self, current_item):
        path = '/' + current_item + '/'
        while len(self.treeview.parent(current_item)) != 0 :
            parent = self.treeview.parent(current_item)
            current_item = parent
            path = '/' + parent + path
        return ('..' + path)

    def search_button_callback(self): # Fonction boutton
        if self.is_in_critical_section:
            return

        self.print_on_window("Searching")
        self.is_in_critical_section = True
        if self.labelText.get() == 'Aucun' :
            self.print_on_window("Nothing selected")
        else :
            self.print_on_window(self.get_path(self.labelText.get()))
            for (index, name, path, NC_list) in self.organism_df.itertuples():
                if path.find(self.get_path(self.labelText.get())) != -1:
                    self.print_on_window("Mais oui c'est cair'")
                    fetch.load_data_from_NC(index, name, path, NC_list)
                    break

        self.window.update()

        self.print_on_window("Research finished")
        self.is_in_critical_section = False

    def reset_button_callback(self):
        if self.is_in_critical_section: return
        progress = ttk.Progressbar(self.window, orient = HORIZONTAL,
            length = 1200, mode = 'determinate')
        progress.pack(side=BOTTOM)
        self.print_on_window("Reset Tree starting")
        self.is_in_critical_section = True
        self.window.update()
        # reset file tree
        self.organism_df = fetch.reset_tree(progress, self.window)
        #reset treeview
        (self.treeview, self.scrollbar) = self.create_tree()
        self.treeview.bind("<<TreeviewSelect>>", self.on_tree_select)
        self.print_on_window("Reset Tree finished")
        progress.destroy()
        self.is_in_critical_section = False

    def on_tree_select(self, event): #on recupere l'organisme dans item
            self.print_on_window("selected items:")
            for item in self.treeview.selection():
                item_text = self.treeview.item(item,"text")
                self.print_on_window(item_text)
                self.labelText.set(item_text)

if __name__ == "__main__":
    App = GUI()