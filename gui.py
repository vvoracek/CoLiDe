from decoot import Decoot, Parameters, Output, BasesToCodons
import data
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from copy import deepcopy as dc
import tkinter as tk
from tkinter import messagebox
from time import time
from PIL import Image, ImageTk
from io import *
from random import shuffle
from save_pdf import save_pdf
import datetime
import warnings
import os
import io
import sys
import platform

class Gui:
    def __init__(self):
        self.root = tk.Tk(className = 'deCOOT')
        self.root.resizable(height = None, width = None)
        W, H = self.root.winfo_screenwidth(), self.root.winfo_screenheight()
        self.SIZE=500
        if(H < 1000):
            if(platform.system() == 'Windows'):
                    self.root.state('zoomed')
            self.SIZE = 300


        plt.figure()

        self.spikedsv = tk.StringVar()
        self.spikedsv.set("2")
        self.lengthsv = tk.StringVar()
        self.removedsv = tk.StringVar()
        self.thresholdsv = tk.StringVar()
        self.reassignsv = tk.StringVar()
        self.inputfilesv = tk.StringVar()
        self.codonsvs = []
        self.text_field_result = None
        self.text_field_stats = None
        self.text_field_string = None
        self.panel = None
        self.output = None
        self.graph_panel = None
        self.button_compute = None

        r = 0
        tk.Radiobutton(self.root, text="Spiked codons     ",variable=self.spikedsv, value="1").\
            grid(row = r, column = 0)
        r += 1
        
        tk.Radiobutton(self.root, text="Degenerate codons ",variable=self.spikedsv, value="2").\
            grid(row = r, column = 0)
        r += 1

        tk.Label(self.root, text="length of AA sequence").grid(row=r, column = 0)

        entry = tk.Entry(self.root, textvariable=self.lengthsv)
        entry.insert(tk.END, '10')
        entry.grid(row = r, column = 1)
        r += 1

        tk.Label(self.root, text="removed tripletes").grid(row=r, column = 0)
        
        entry = tk.Entry(self.root, textvariable=self.removedsv)
        entry.insert(tk.END, '')
        entry.grid(row = r, column = 1)
        r += 1

        tk.Label(self.root, text="maximum rate").grid(row=r, column = 0)
        
        entry = tk.Entry(self.root, textvariable=self.thresholdsv)
        entry.insert(tk.END, '0.9')
        entry.grid(row = r, column = 1)
        r += 1

        tk.Label(self.root, text="model distribution").grid(row=r, column = 0)
        self.listbox = tk.Listbox(self.root, height = 7, selectmode = 'SINGLE')
        for option in data.options:
            self.listbox.insert(tk.END, option)
        self.listbox.selection_set(0)

        self.listbox.insert(tk.END, "        none")
        self.listbox.grid(row = r, column = 1)
        r += 1

        self.listbox.configure(exportselection=False)

        tk.Label(self.root, text="reassign codons").grid(row=r, column = 0)
        
        entry = tk.Entry(self.root, textvariable=self.reassignsv)
        entry.grid(row = r, column = 1)
        r += 1
        vec2fit = [.1208,.0,.0,.0734,.0,.0,.0,.0,.0,.0905,.0679,.118,.122,.157,.0,.0,.0388,.0394,.1722,.0,.0] 
        for i in range(len(data.names)):
            self.codonsvs.append(tk.StringVar())
            self.codonsvs[-1].set(vec2fit[i])
            s = ""
            if(type(data.colors[i]) is str):
                color = data.colors[i]
            else:
                color = "#"
                c1 = data.colors[i][0]
                c2 = data.colors[i][1]
                c3 = data.colors[i][2]
                color = '#' + ('00' + hex(int(255*c1))[2:])[-2:] + ('00' + hex(int(255*c2))[2:])[-2:] + ('00' + hex(int(255*c3))[2:])[-2:]
                
            tk.Label(self.root, text=data.names[i] + s, bg = color, fg = "#ffffff", font=('Courier', 12)).grid(row=r, column = 0)
            tk.Entry(self.root, textvariable=self.codonsvs[-1]).grid(row = r, column = 1)
            r += 1
       


        img2 = Image.new("RGBA", (self.SIZE, self.SIZE))
        img2 = ImageTk.PhotoImage(img2)

        self.panel = tk.Label(self.root, image=img2, width = self.SIZE, height = self.SIZE)
        self.panel.grid(row = 16, column = 2, rowspan = 100)
   

        self.text_field_result = tk.Text(self.root, height=25, width=39)
        s = "     expected   reached  difference\n"
        for i, j, k in zip([0]*21, [0]*21, range(21)):
            if(len(data.names[k]) == 3):
                s += data.names[k] + "    0         0         0\n"
            else:
                s += data.names[k] + "   0         0         0\n"

        self.text_field_result.insert(tk.END, s)
        self.text_field_result.grid(row = 1, rowspan = 17, column = 4)

        img2 = Image.new("RGBA", (self.SIZE//10*7, self.SIZE))
        img2 = ImageTk.PhotoImage(img2)
        self.graph_panel = tk.Label(self.root, image=img2, width=self.SIZE, height = (self.SIZE*7)//10)
        self.graph_panel.grid(row = 1, column = 2, rowspan = 15)

        self.text_field_stats = tk.Text(self.root, height=8, width=39)
        s = "mean error of a single protein\nunknown\n"
        s +="variance of error of a single protein\nuknown\n"
        s +="mean GC content of the DNA template\nuknown\n"
        s +="mean mass of a protein\nuknown"

        self.text_field_stats.insert(tk.END, s)
        self.text_field_stats.grid(row = 19, rowspan = 8, column = 4)

        self.text_field_string = tk.Text(self.root, height = 1, width = 39)
        self.text_field_string.insert(tk.END, "output string")
        self.text_field_string.grid(row = 28, column = 4)

        tk.Label(self.root, text="input from file").grid(row=r, column = 0)
        
        entry = tk.Entry(self.root, textvariable=self.inputfilesv)
        entry.grid(row = r, column = 1)
        r += 1

        tk.Button(self.root, text='clear', command = lambda:[self.codonsvs[i].set(0) for i in range(len(vec2fit))]+[self.reassignsv.set("")]).grid(row=r, column=1)
        self.button_compute = tk.Button(self.root, text="compute", command=self.callback_compute)
        self.button_compute.grid(row = r+1, column = 1)


        tk.Button(self.root, text="permute codons", command=self.permute_codons).grid(row = 50, column=1)
        tk.Button(self.root, text="export to pdf & save imgs", command=self.callback_export).grid(row = 50, column=4)
        tk.Label(self.root).grid(row = r+4, column = 5)  
        self.root.protocol("WM_DELETE_WINDOW", sys.exit)
        for row_num in range(self.root.grid_size()[1]):
              self.root.rowconfigure(row_num, minsize=0, weight=1)
        for col_num in range(self.root.grid_size()[0]):
              self.root.columnconfigure(col_num, minsize=0, weight=1)


        self.root.mainloop()

    def callback_export(self):
        now = datetime.datetime.now()
        name = "decoot_" + str(now).replace(" ", "_").replace(".", "-").replace(":", ";")
        self.output.img.save("sequence_" + name + ".png")
        self.output.graph_error.save("graph_" + name + ".png")
        save_pdf(self.output, name)
    


    def callback_compute(self):
        self.button_compute.config(state='disabled')
        self.output = None
        self.clear_img()
        self.text_field_stats.delete('1.0', tk.END)
        self.text_field_result.delete('1.0', tk.END)
        self.text_field_string.delete('1.0', tk.END)

        vec2fit = 21*[0]
        if(self.inputfilesv.get().split('.')[-1] == 'txt'):
            [self.codonsvs[i].set(0) for i in range(len(vec2fit))]
            dict_of_params = {}
            self.spikedsv.set('2')
            with open(self.inputfilesv.get(), 'r' ) as f:
                for line in f:
                    if('spiked' in line):
                        self.spikedsv.set('1')
                    if('=' in line):
                        if('\n' in line):
                            line = line.split('\n')[0]
                        line = line.split('=')
                        dict_of_params[line[0]]=line[1]

            for i in range(21):
                if(data.names[i].split()[0] in dict_of_params):
                    self.codonsvs[i].set(dict_of_params[data.names[i].split()[0]])


            
            if('model_distribution' in dict_of_params):
                model_distribution = None
                mod = dict_of_params['model_distribution']
                idx = 0
                self.listbox.selection_clear(0, tk.END)

                for opt in data.options:
                    if(mod in opt.split()[0]):
                        self.listbox.selection_set(idx)
                        break
                    idx += 1
                self.listbox.selection_set(idx)

            if('maximum_rate' in dict_of_params):
                threshold = dict_of_params['maximum_rate']
                self.thresholdsv.set(str(threshold))
            else:
                threshold = float(self.thresholdsv.get())

            if('length' in dict_of_params):
                length = dict_of_params['length']
                self.lengthsv.set(str(length))
            else:
                length = int(self.lengthsv.get())

            if('removed_triplets' in dict_of_params):
                self.removedsv.set(dict_of_params['removed_triplets'])

            if('reassigned_codons' in dict_of_params):
                self.reassignsv.set(dict_of_params['reassigned_codons'])

                
        for i, j in zip(range(21), self.codonsvs):
            vec2fit[i] = float(j.get())

        reas = [i.split() for i in self.reassignsv.get().split(',')]
        if(not reas[-1]):
            reas = reas[:-1]
        comb = sum(vec2fit) + sum([float(i[1])for i in reas])
        
        tmp = ""
        for i in reas:
            tmp += i[0] + " " + str(round(float(i[1])/comb, 3)) + ", "
        self.reassignsv.set(tmp[:-2])

        reas = [i.split() for i in self.reassignsv.get().split(',')]
        if(not reas[-1]):
            reas = reas[:-1]

        if(vec2fit is not None and comb > 0 and min(vec2fit)>=0):
            vec2fit= [i/comb for i in vec2fit]
        else:
            messagebox.showerror("Error", "Invalid codon rates.")
            self.button_compute.config(state='normal')
            return

        for i in range(len(vec2fit)):
            self.codonsvs[i].set(round(vec2fit[i], 3))


        b2c = BasesToCodons()
        for r in reas:
            for codon in b2c.codons:
                if(r[0] in codon):
                    codon.remove(r[0])
                    b2c.codons.append([r[0]])
                    vec2fit.append(float(r[1]))
                    break

        rem = self.removedsv.get().upper().split()
        transposition_table = {'T' : 0, 
                            'C' : 1, 
                            'A' : 2, 
                            'G' : 3}
        removed_triplets = []
        s = ""
        rmv = []
        for i in rem:
            if(len(i) == 3 and \
            i[0] in transposition_table and \
            i[1] in transposition_table and \
            i[2] in transposition_table ):
                tmp = [[0,0,0,0], [0,0,0,0], [0,0,0,0]]
                for j in range(3):
                    tmp[j][transposition_table[i[j]]] = 1
                rmv.append(i)
                removed_triplets.append(tmp)

        spiked_codons = self.spikedsv.get() == "1"
        try:
            model_distribution = data.options[self.listbox.get(self.listbox.curselection()[0])]
        except:
            model_distribution = None

        length = int(self.lengthsv.get())
        threshold = float(self.thresholdsv.get())
        if(length < 1):
            messagebox.showerror("Error", "Invalid length of sequence.")
            self.button_compute.config(state='normal')
            return
        elif(threshold > 1 or threshold < 0):
            messagebox.showerror("Error", "Maximum rate should be in [0,1].")
            self.button_compute.config(state='normal')
            return
        self.root.update()
        parameters = Parameters(model_distribution, threshold, spiked_codons, removed_triplets, vec2fit, length, b2c)
        decoot = Decoot(parameters)
        output = decoot()
        self.output = output
        self.output.removed_triplets = rmv
        if(not self.output):
            messagebox.showerror("Error", "Could not find a feasible solution.")
            self.button_compute.config(state='normal')
            return
        self.output.make_imgs()
        
        img = self.output.graph_error.resize((self.SIZE, (self.SIZE*7)//10), Image.ANTIALIAS)
        img = ImageTk.PhotoImage(img)
        self.graph_panel.configure(image = img)
        self.graph_panel.image = img

        s, v = self.output.img.size
        w = self.SIZE
        h = int(w*v/s)
        img = self.output.img.resize((w, h), Image.ANTIALIAS)
        img = ImageTk.PhotoImage(img)
        self.panel.configure(image = img)
        self.panel.image=img

        self.text_field_string.delete('1.0', tk.END)
        self.text_field_string.insert(tk.END, output.output_string)

        self.text_field_stats.delete('1.0', tk.END)
        s = "mean error of a single protein\n"
        s += str(self.output.mean) + "\n"
        s +="variance of error of a single protein\n"
        s += str(self.output.var) + "\n"
        s +="mean GC content of a protein\n"
        if(output.vec2fit[20] != 0):
            self.output.gc = 'mixed species'
        s += str(self.output.gc) + "%\n"
        s +="mean mass of a protein\n"
        s += str(self.output.weight)
        self.text_field_stats.insert(tk.END, s)

        self.text_field_result.delete('1.0', tk.END)
        s = "       expected   reached  difference\n"
        for i, j, k in zip(vec2fit, self.output.reached_distribution, data.names[:-1] + [" * "] + sum(self.output.parameters.b2c.codons[21:], [])):
            s += k + "    "
            if(i):
                s += "%.3f    " %(i)
            else:
                s += "0        "
            if(j):
                s += "%.3f      " %(j/sum(self.output.reached_distribution))
            else:
                s += "0          "
            if(abs(i-j/sum(self.output.reached_distribution)) > 10**-5):
                if(i-j/sum(self.output.reached_distribution) < 0):
                    s += " "
                s += " %.3f\n" %(round(-(i-j/sum(self.output.reached_distribution)), 5))
            else:
                s += "  0\n"
    
        self.text_field_result.insert(tk.END, s)
        self.button_compute.config(state='normal')



    def clear_img(self):
        img = Image.new("RGBA", (self.SIZE, self.SIZE))
        img = ImageTk.PhotoImage(img)
        self.panel.configure(image = img)
        self.panel.image=img

        img = Image.new("RGBA", ((self.SIZE*6)//10,(self.SIZE*2)//5))
        img = ImageTk.PhotoImage(img)
        self.graph_panel.configure(image = img)
        self.graph_panel.image = img

    def permute_codons(self):
        if(not self.output):
            messagebox.showerror("Error", "No sequence to permute.")
            return
        self.output.shuffle()
        s, v = self.output.img.size
        w = self.SIZE
        h = int(w*v/s)
        img = self.output.img.resize((w, h), Image.ANTIALIAS)
        img = ImageTk.PhotoImage(img)
        self.panel.configure(image = img)
        self.panel.image=img

if(__name__ == '__main__'):
    Gui()
