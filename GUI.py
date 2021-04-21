from tkinter import *
from tkinter import ttk

from run import *
import os
import time


# Creating window
window = Tk()
window.geometry("1200x800")
window.title("Michel on t'adore")



# Frame 1 : Liste des fichiers
list = Frame(window)
list.place(x=0, y=0, anchor="nw", width=400, height=800)

# Creating treeview window
scrollbar = Scrollbar(list)
scrollbar.pack( side = RIGHT, fill = Y )

treeview = ttk.Treeview(list)
treeview.heading('#0', text='Arborescence des fichiers')

treeview.configure(yscrollcommand=scrollbar.set)
treeview.pack(fill="both", expand=True)
scrollbar.configure(command=treeview.yview)


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

Label(Info_lab, text="Organisme choisi : ").grid(row = 0, column = 0, sticky = W, ipadx = 100, ipady = 30)

labelText = StringVar()
labelText.set('Aucun')
depositLabel = Label(Info_lab, textvariable=labelText)
depositLabel.grid(row = 0, column = 1, sticky = W)

Label(Info_lab, justify = LEFT, text="Région fonctionnelle choisie : ").grid(row = 4, column = 0, sticky = W, ipadx = 100, ipady = 30)



# Frame 3 : Logs
Logs = Frame(window, background="#b22222")
Logs.place(x=400, y=400, anchor="nw", width=800, height=400)

Logs_lab = LabelFrame(Logs, text="Logs")
Logs_lab.pack(fill="both", expand="yes")


log_text = Text(Logs_lab, height='400', width='800')
scroll = Scrollbar(Logs_lab, command = log_text.yview)
log_text.configure(yscrollcommand=scroll.set)

log_text.pack(side=LEFT)
scroll.pack(side=RIGHT, fill=Y)

def print_on_window(t): #affiche t dans les logs
    time_string = time.strftime('%H:%M:%S')
    log_text.insert(INSERT, time_string + ' : ' + t + "\n")
    log_text.yview(END)




# Menu des regions

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


variable = StringVar(Info_lab)
variable.set(OptionList[0])

menu = OptionMenu(Info_lab, variable, *OptionList)
menu["borderwidth"]=1
menu.grid(row = 4, column = 1, sticky = W)


def callback(*args): # fonction pour executer du code pour le menu, a changer
    print_on_window("The selected item is "+variable.get())

variable.trace("w", callback)

def button_callback(): # Fonction boutton
   print_on_window("Run...")

run = Button(Info_lab, text ="Run", command = button_callback, relief = RIDGE, borderwidth=1)
run.grid(row = 5, sticky = 'se', column = 2, ipadx = 20, pady = 30, padx = 30)

# Fonctionnalités

def on_tree_select(event): #on recupere l'organisme dans item
        print_on_window("selected items:")
        for item in treeview.selection():
            item_text = treeview.item(item,"text")
            print_on_window(item_text)
            labelText.set(item_text)


treeview.bind("<<TreeviewSelect>>", on_tree_select)



window.mainloop()
