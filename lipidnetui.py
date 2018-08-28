# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import tkinter as tk
import time
from tkinter import ttk
from Bio import Entrez
from bs4 import BeautifulSoup as soup
Entrez.email = "md.nur.hakim.rosli@gmail.com"
    
class LipidUI(tk.Tk):

    def __init__(self, *args, **kwargs):
        base = tk.Tk.__init__(self, *args, **kwargs)
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0,weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}
        for pg in (HomePage, SearchPage, FetchPage):
            print(pg)
            page = str(pg).replace("<class '__main__.", "")
            page = page.replace("'>","")
            frame = pg(container, self)
            self.frames[page] = frame
            frame.grid(row=0, column=0, sticky="nsew")
                
        self.show_frame('HomePage')

    def show_frame(self, cont):
        cont = str(cont).replace("<class '__main__.", "")
        Cont = cont.replace("'>","")
        frame = self.frames[Cont]
        frame.tkraise()

            
class HomePage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="WELCOME TO LIPIDNET", font='LARGE_FONT')
        label.pack(padx=10, pady=10)
        
        button1 = tk.Button(self, text="SEARCH",  bg='gold', fg='yellow',
                           command=lambda: controller.show_frame('SearchPage'))
        button1.pack(padx=10, pady=10, side='left')
        button2 = tk.Button(self, text="FETCH",  bg='gold', fg='yellow',
                           command=lambda: controller.show_frame('FetchPage'))
        button2.pack(padx=10, pady=10, side='right')
        
class SearchPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        entry = tk.Entry(self)
        title = tk.Label(self, text="LIPIDNET SEARCH PAGE", font='LARGE_FONT')
        title.pack(padx=10, pady=10)
        label = tk.Label(self, text="Please enter search term below:", font='LARGE_FONT')
        label.pack(padx=10, pady=10)
        entry = tk.Entry(self)
        entry.pack(padx=10, pady=10)
        self.term = entry.get()
        self.progress = ttk.Progressbar(self, orient="horizontal",
                                        length=200, mode="determinate")
        searchbutton = tk.Button(self, text="Obtain records", command=self.start)
        searchbutton.pack(side="top", padx=10, pady=10)
        
        self.progress = ttk.Progressbar(self, orient="horizontal",
                                        length=200, mode="determinate")
        self.progress.pack(side="top")
        
        button1 = tk.Button(self, text="HOME",  bg='gold', fg='yellow',
                           command=lambda: controller.show_frame('HomePage'))
        button1.pack(padx=10, pady=10, side='left')
        button2 = tk.Button(self, text="FETCH",  bg='gold', fg='yellow',
                           command=lambda: controller.show_frame('FetchPage'))
        button2.pack(padx=10, pady=10, side='right')
        
    def start(self):
        self.progress["value"] = 0
        self.progress["maximum"] = 100
        self.start = time.time()
        
        label = tk.Label(self, text="Beginning Search", font='LARGE_FONT')
        label.pack(padx=10, pady=10)

        s = Entrez.esearch(db="pubmed", term=self.term, retmode="xml").read()
        print(s)
        idlist = [id.text for id in s.find_all("Id")]

        self.update()
        self.timetaken = self.end - self.start
        timeresult = "Search completed in {} seconds".format(self.timetaken)
        full_count = "{} records obtained from a total of {}.".format(s.RetMax.text, s.Count.text)
        result = tk.Label(self, text=timeresult)
        count = tk.Label(self, text=full_count)
        count.pack(padx=10, pady=10)
        result.pack(padx=10, pady=10)
    
    def update(self):
        self.progress["value"] = 100
        self.end = time.time()
        
class FetchPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent, bg='yellow')
        label = tk.Label(self, text="LIPIDNET FETCH PAGE", font='LARGE_FONT', bg='yellow', fg='gold')
        label.pack(padx=20, pady=20)
        self.progress = ttk.Progressbar(self, orient="horizontal",
                                        length=200, mode="determinate")
        button1 = tk.Button(self, text="HOME", bg='gold', fg='yellow',
                           command=lambda: controller.show_frame(HomePage))
        button1.pack(padx=10, pady=10, side='left')
        button2 = tk.Button(self, text="SEARCH", bg='gold', fg='yellow',
                           command=lambda: controller.show_frame(SearchPage))
        button2.pack(padx=10, pady=10, side='right')
        self.progress.pack()
        
    def start(self):
        self.progress["value"] = 0
        self.progress["maximum"] = 100
        self.start = time.time()
        idlist = Entrez.read(Entrez.esearch(db="pubmed", term=self.term))
        label = tk.Label(self, text="Beginning Search", font=LARGE_FONT)
        label.pack(padx=10, pady=10)
        self.update()
        self.timetaken = self.end - self.start
        timeresult = "Search completed in {} seconds".format(self.timetaken)
        result = tk.Label(self, text=timeresult)
        result.pack(padx=10, pady=10)
    
    def update(self):
        self.progress["value"] = 100
        self.end = time.time()
             
        
        
    

app = LipidUI()
app.mainloop()