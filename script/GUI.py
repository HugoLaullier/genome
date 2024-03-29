#!/usr/bin/env python3
from tkinter import *
from tkinter import ttk
from ttkwidgets import CheckboxTreeview

import os
import time
import ctypes

import fetch as fetch

VERBOSE = False

try:
    ctypes.windll.shcore.SetProcessDpiAwareness(True) # améliore la netteté de l'app
except:
    pass # doesn't work on linux

def launch_welcome() :
    root = Tk()
    root.title("Genbank - Welcome")

    Label(root, text="Création de l'arborescence en cours...").pack(side = LEFT, anchor=  W, padx = 15, pady = 15)
    progressbar = ttk.Progressbar(root, orient = HORIZONTAL, length = 100, mode = 'determinate')
    progressbar.pack( padx = 15, pady = 15)

    root.after(200,  lambda : fetch.reset_tree(root, progressbar))
    root.mainloop()

class GUI:
    def __init__(self):
        # attribute
        launch_welcome()
        self.organism_df = fetch.load_df_from_pickle()

        # Creating window
        self.window = Tk()
        self.window.geometry("1200x800")
        self.window.title("GenBank Application")
        
        self.tree_array = []
        self.org_selected = []
        # Frame 1 : Liste des fichiers
        (self.treeview, self.scrollbar) = self.create_tree()
        #self.treeview.bind("<<TreeviewSelect>>", self.on_tree_select)
        self.treeview.bind("<Button-1>", self.box_click, True)
        # Frame 2 : Informations
        (self.Info, self.Info_lab1, self.labelText) = self.create_info()

        # Frame 3 : Logs
        (self.Logs, self.Logs_lab, self.log_text, self.log_scroll) = self.create_log()

        # Regions menu
        (self.selected_region, self.menu_menu, self.run_search, self.OptionList_region) = self.create_region_menu()

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
        name_nodes.sort()
        if name_nodes == [] : pass
        for name in name_nodes :
            name_path = node_path + name
            if os.path.isdir(name_path):
                name = self.tree_exist(tree, name)
                try :
                    tree.insert(node, '1000000',iid=name, text=name)
                    self.tree_array.append(name)
                    name_path += '/'
                    self.create_node(tree, name_path, name)
                except:
                    if VERBOSE:
                        print(name_path)
                    pass

    def create_tree(self):
        list = Frame(self.window)
        list.place(x=0, y=0, anchor="nw", width=400, height=800)

        # Creating treeview window
        scrollbar = Scrollbar(list)
        scrollbar.pack( side = RIGHT, fill = Y )

        treeview = CheckboxTreeview(list)
        treeview.heading('#0', text='Arborescence des fichiers')
        treeview.configure(yscrollcommand=scrollbar.set)
        treeview.pack(fill="both", expand=True)
        
        scrollbar.configure(command=treeview.yview)

        root_path = "../Results/"
        treeview.insert('', '0', text='Results', iid='Results')
        
        self.tree_array.append('Results')
        self.create_node(treeview, root_path, 'Results')

        return treeview, scrollbar

    def box_click(self, event):
        """ check or uncheck box when clicked """
        x, y, widget = event.x, event.y, event.widget
        elem = widget.identify("element", x, y)
        if "image" in elem:
            # a box was clicked
            item = self.treeview.identify_row(y)
            tags = self.treeview.item(item, "tags")
            if ("unchecked" in tags) or ("tristate" in tags):
                try : 
                    self.org_selected.remove(item)
                except : 
                    pass
            else:
                self.org_selected.append(item)

    # INFO CREATION METHODS
    def create_info(self):
        Info = Frame(self.window)
        Info.place(x=400, y=0, anchor="nw", width=800, height=150)

        Info_lab1 = LabelFrame(Info, text="Selectionner une région fonctionnelle", padx=20, pady=20)
        Info_lab1.pack(fill="both", expand="yes")

        labelText = StringVar(value='Aucun')

        Label(Info_lab1, justify = LEFT, text="Région fonctionnelle choisie : ").grid(row = 1, column = 0, sticky = W, ipadx = 50, ipady = 10)
        return Info, Info_lab1, labelText

    # LOG CREATION METHODS
    def create_log(self):
        Logs = Frame(self.window, background="#b22222")
        Logs.place(x=400, y=150, anchor="nw", width=802, height=650)

        Logs_lab = LabelFrame(Logs, text="Logs")
        Logs_lab.pack(fill="both", expand="yes")


        log_text = Text(Logs_lab, height=400, wrap=WORD, width=800, padx=5, pady=5)
        log_scroll = Scrollbar(Logs_lab, command = log_text.yview)
        log_text.configure(yscrollcommand=log_scroll.set)
        log_text.bindtags((str(log_text), str(Logs_lab), "all"))
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

        variable = StringVar(self.Info_lab1)
        variable.set(OptionList[0])

        menu = OptionMenu(self.Info_lab1, variable, *OptionList)
        menu["borderwidth"]=1
        menu.grid(row = 1, column = 1, sticky = W, pady = 15)

        variable.trace("w", self.callback)

        run_search = Button(self.Info_lab1, text ="Search", command = self.search_button_callback, borderwidth=1)
        run_search.grid(row = 1, column = 2, ipadx = 20, pady = 15, padx = 50)
        return variable, menu, run_search, OptionList

    # FEATURES METHODS

    def update_tree_tags(self):
        for node in self.tree_array:
            current_path = self.get_path(node).replace("Other1", "Other").replace("unclassified1", "unclassified").replace("uncultured_bacterium1", "uncultured_bacterium")
            is_leaf = False
            for (index, name, path, NC_list) in self.organism_df.itertuples():
                path_full = path + name.replace(" ", "_").replace("[", "_").replace("]", "_").replace(":", "_") + '/'
                if path_full == current_path:
                    is_leaf = True
                    if os.listdir(path_full) == []:
                        self.treeview.item(node)
                    else:
                        self.treeview.item(node)
                    break
            if is_leaf:
                continue

        self.print_on_window("Tree updated")

    def print_on_window(self, t): #affiche t dans les logs
        time_string = time.strftime('%H:%M:%S')
        self.log_text.insert(INSERT, time_string + ' : ' + str(t) + "\n")
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

        self.print_on_window("Searching")
        if self.org_selected == [] :
            self.print_on_window("No organism selected")
            self.window.update()
            return
        elif self.selected_region.get() == 'Aucun':
            self.print_on_window("No functional region selected")
            self.window.update()
            return
        else :
            org_done = []
            for org in self.org_selected :
                current_path = self.get_path(org).replace("Other1", "Other").replace("unclassified1", "unclassified").replace("uncultured_bacterium1", "uncultured_bacterium")
                c = 0
                nb_region_found = 0
                for (index, name, path, NC_list) in self.organism_df.itertuples():
                    path_full = path + name.replace(" ", "_").replace("[", "_").replace("]", "_").replace(":", "_") + '/'
                    nb_new_region_found = -1
                    if path_full == current_path and name not in org_done:
                        if VERBOSE:
                            print(name)
                        c += 1
                        nb_region_found_1 = 0
                        nb_region_already_downloaded_1 = 0
                        if(self.selected_region.get() == 'CDS'):
                            nb_region_found_1, nb_region_already_downloaded_1 = fetch.load_data_from_NC(index, name, path, NC_list, 'intron')
                        self.window.update()
                        nb_region_found_2, nb_region_already_downloaded_2 = fetch.load_data_from_NC(index, name, path, NC_list, self.selected_region.get())
                        nb_new_region_found = nb_region_found_1 + nb_region_found_2
                        nb_new_region_already_downloaded = nb_region_already_downloaded_1 + nb_region_already_downloaded_2
                        org_done.append(name)
                        if nb_new_region_found == 0 and nb_new_region_already_downloaded == 0:
                            self.print_on_window("Selected functional region [" + self.selected_region.get() + "] not found for organism [" + name + "]")
                            c -= 1
                        elif nb_new_region_found == 0 and nb_new_region_already_downloaded != 0:
                            self.print_on_window("[" + self.selected_region.get() + "] for [" + name + "] already downloaded ")
                        else:
                            self.print_on_window("[" + self.selected_region.get() + "] for [" + name + "] downloaded ")
                        self.window.update()
                        nb_region_found += nb_new_region_found
                        break
                    elif current_path in path_full and name not in org_done:
                        c += 1
                        nb_region_found_1 = 0
                        nb_region_already_downloaded_1 = 0
                        self.window.update()
                        if(self.selected_region.get() == 'CDS'):
                            nb_region_found_1, nb_region_already_downloaded_1 = fetch.load_data_from_NC(index, name, path, NC_list, 'intron')
                        nb_region_found_2, nb_region_already_downloaded_2 = fetch.load_data_from_NC(index, name, path, NC_list, self.selected_region.get())

                        nb_new_region_found = nb_region_found_1 + nb_region_found_2
                        nb_new_region_already_downloaded = nb_region_already_downloaded_1 + nb_region_already_downloaded_2

                        org_done.append(name)
                    if nb_new_region_found != -1:
                        if VERBOSE:
                            print(name)

                        if nb_new_region_found == 0 and nb_new_region_already_downloaded == 0:
                            self.print_on_window("Selected functional region [" + self.selected_region.get() + "] not found for organism [" + name + "]")
                            c -= 1
                        elif nb_new_region_found == 0 and nb_new_region_already_downloaded != 0:
                            self.print_on_window("[" + self.selected_region.get() + "] for [" + name + "] already downloaded ")
                        else:
                            self.print_on_window("[" + self.selected_region.get() + "] for [" + name + "] downloaded ")
                        self.window.update()
                        nb_region_found += nb_new_region_found
                # if nb_region_found == 0:
                #     print("hello")
                self.window.update()


        self.print_on_window(str(c) + " items downloaded")
        self.window.update()
        self.print_on_window("--------------------Research finished--------------------")


    def on_tree_select(self, event): #on recupere l'organisme dans item, A GARDER ?????
        for item in self.treeview.selection():
            item_text = self.treeview.item(item,"text")
            self.labelText.set(item_text)

if __name__ == "__main__":
    App = GUI()