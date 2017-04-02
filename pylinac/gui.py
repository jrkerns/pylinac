from tkinter import *
from tkinter.ttk import *
from tkinter import filedialog
import os.path as osp
import os

from pylinac import PicketFence, vmat


class PylinacGUI(Frame):

    def __init__(self, master=None):
        super().__init__(master)
        self.pack()
        self.notebook = Notebook(self)
        self.init_pf()
        self.init_vmat()
        for child in self.winfo_children():
            child.grid_configure(padx=10, pady=10)

    def init_vmat(self):

        def load_open():
            f = filedialog.askopenfilename()
            self.vmat_openimg.set(f)

        def load_dmlc():
            f = filedialog.askopenfilename()
            self.vmat_dmlcimg.set(f)

        def analyze_vmat():
            v = vmat.VMAT(images=(self.vmat_openimg.get(), self.vmat_dmlcimg.get()), delivery_types=('open', 'dmlc'))
            v.analyze(test=self.vmat_test.get(), tolerance=self.vmat_tol.get(), x_offset=self.vmat_xoff.get())
            fname = osp.join(self.vmat_dmlcimg.get().replace('.dcm', '.pdf'))
            v.publish_pdf(fname)
            self.vmat_pdf.set(fname)
            os.startfile(fname)

        self.vmat_tab = Frame(self.notebook)
        self.vmat_openimg = StringVar()
        self.vmat_dmlcimg = StringVar()
        self.vmat_test = StringVar(value=vmat.DRGS)
        self.vmat_tol = DoubleVar(value=1.5)
        self.vmat_xoff = IntVar(value=0)
        self.vmat_pdf = StringVar()
        Button(self.vmat_tab, text='Load Open Image...', command=load_open).grid(column=1, row=1)
        Button(self.vmat_tab, text='Load DMLC Image...', command=load_dmlc).grid(column=1, row=3)
        Label(self.vmat_tab, textvariable=self.vmat_openimg).grid(column=1, row=2)
        Label(self.vmat_tab, textvariable=self.vmat_dmlcimg).grid(column=1, row=4)
        Combobox(self.vmat_tab, values=(vmat.DRGS, vmat.DRMLC), textvariable=self.vmat_test).grid(column=1, row=5)
        Label(self.vmat_tab, text='Tolerance (%):').grid(column=1, row=6)
        Entry(self.vmat_tab, width=7, textvariable=self.vmat_tol).grid(column=2, row=6)
        Label(self.vmat_tab, text='X-offset (px):').grid(column=1, row=7)
        Entry(self.vmat_tab, width=7, textvariable=self.vmat_xoff).grid(column=2, row=7)
        Button(self.vmat_tab, text='Analyze', command=analyze_vmat).grid(column=1, row=8)
        Label(self.vmat_tab,
              text='Analysis will analyze the file(s) according to the settings, \nsave a PDF in the same directory as the original file location and then open it.').grid(
            column=1, row=9)
        Label(self.vmat_tab, text='Save file:').grid(column=1, row=10)
        Label(self.vmat_tab, textvariable=self.vmat_pdf).grid(column=1, row=11)
        self.notebook.add(self.vmat_tab, text='VMAT')
        for child in self.vmat_tab.winfo_children():
            child.grid_configure(padx=10, pady=5)

    def init_pf(self):

        def load_file():
            f = filedialog.askopenfilename()
            self.pf_file.set(f)

        def analyze_pf():
            if self.pf_filter.get():
                pf = PicketFence(self.pf_file.get(), filter=3)
            else:
                pf = PicketFence(self.pf_file.get())
            atol = self.pf_atol.get() if self.pf_atol.get() == 0 else None
            pickets = self.pf_pickets.get() if self.pf_pickets.get() == 0 else None
            hd = self.pf_hdmlc.get()
            pf.analyze(tolerance=self.pf_tol.get(),
                       action_tolerance=atol,
                       hdmlc=hd,
                       num_pickets=pickets,
                       )
            fname = osp.join(self.pf_file.get().replace('.dcm', '.pdf'))
            pf.publish_pdf(fname)
            os.startfile(fname)

        self.pf_tab = Frame(self.notebook)
        self.pf_filter = BooleanVar(value=False)
        self.pf_file = StringVar()
        self.pf_tol = DoubleVar(value=0.5)
        self.pf_atol = DoubleVar(value=0.25)
        self.pf_pickets = IntVar(value=10)
        self.pf_hdmlc = BooleanVar(value=False)
        Checkbutton(self.pf_tab, text='Filter upon loading', variable=self.pf_filter).grid(column=1, row=1)
        Button(self.pf_tab, text='Load File...', command=load_file).grid(column=1, row=2)
        Label(self.pf_tab, text='File:').grid(column=1, row=3)
        Label(self.pf_tab, textvariable=self.pf_file).grid(column=2, row=3)
        Label(self.pf_tab, text='Tolerance (mm):').grid(column=1, row=4)
        Entry(self.pf_tab, width=7, textvariable=self.pf_tol).grid(column=2, row=4)
        Label(self.pf_tab, text='(Optional) Action Tolerance (mm):').grid(column=1, row=5)
        Entry(self.pf_tab, width=7, textvariable=self.pf_atol).grid(column=2, row=5)
        Label(self.pf_tab, text='(Optional) Number of pickets:').grid(column=1, row=6)
        Entry(self.pf_tab, width=7, textvariable=self.pf_pickets).grid(column=2, row=6)
        Checkbutton(self.pf_tab, text='HDMLC?:', variable=self.pf_hdmlc).grid(column=1, row=7)
        Button(self.pf_tab, text='Analyze', command=analyze_pf).grid(column=1, row=8)
        Label(self.pf_tab, text='Analysis will analyze the file according to the settings, \nsave a PDF in the same directory as the original file location and then open it.').grid(column=1, row=9)
        self.notebook.add(self.pf_tab, text='Picket Fence')
        for child in self.pf_tab.winfo_children():
            child.grid_configure(padx=10, pady=5)


def gui():
    root = Tk()
    root.title('Pylinac GUI')
    app = PylinacGUI(master=root)
    app.mainloop()
