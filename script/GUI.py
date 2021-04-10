from tkinter import *
from tkinter import ttk

from run import *
import os


# Creating window
window = Tk()
window.geometry("1200x800")
window.title("Michel on t'adore")



# Frame 1 : Liste des fichiers
list = Frame(window)
list.place(x=0, y=0, anchor="nw", width=400, height=1200)

# Creating treeview window
treeview = ttk.Treeview(list)
treeview.heading('#0', text='Arborescence des fichiers')

treeview.pack(fill="both", expand=True)


def tree_exist(tree, name) :
    if tree.exists(name) :
        name += '1'
        tree_exist(tree, name)
    return name



def create_node(tree, node_path, node) :
    name_nodes = os.listdir(node_path)
    if name_nodes == [] : pass
    for name in name_nodes :
        name_path = node_path + name
        if os.path.isdir(name_path) :
            name = tree_exist(tree, name)
            name_path += '/'
            tree.insert(node, '1000000',iid=name, text=name)
            create_node(tree, name_path, name)


root_path = "../Results/"
treeview.insert('', '0', text='Result', iid='Result')
create_node(treeview, root_path, 'Result')


# Frame 2 : Informations
Info = Frame(window)
Info.place(x=400, y=0, anchor="nw", width=800, height=600)

Info_lab = LabelFrame(Info, text="Informations", padx=20, pady=20)
Info_lab.pack(fill="both", expand="yes")

Label(Info_lab, justify = LEFT, text="Organisme choisi : ").grid(row = 0, column = 0)

labelText = StringVar()
labelText.set('Aucun')
depositLabel = Label(Info_lab, textvariable=labelText)
depositLabel.grid(row = 0, column = 1)

Label(Info_lab, justify = LEFT, text="Région fonctionnelle choisie : ").grid(row = 4, column = 0)

# Frame 3 : Logs
Logs = Frame(window, background="#b22222")
Logs.place(x=400, y=400, anchor="nw", width=800, height=400)

Logs_lab = LabelFrame(Logs, text="Logs", padx=20, pady=20)
Logs_lab.pack(fill="both", expand="yes")




# Fonctionnalités

def on_tree_select(event):
        print("selected items:")
        for item in treeview.selection():
            item_text = treeview.item(item,"text")
            print(item_text)
            labelText.set(item_text)


treeview.bind("<<TreeviewSelect>>", on_tree_select)



window.mainloop()
