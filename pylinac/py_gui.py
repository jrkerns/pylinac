import itertools
from tkinter import Label, Button, StringVar, BooleanVar, DoubleVar, IntVar, Checkbutton, Tk
from tkinter.ttk import Frame, Notebook, Combobox, Entry
from tkinter import filedialog
from tkinter import messagebox
import os.path as osp
import webbrowser

from . import picketfence, vmat, ct, log_analyzer, starshot, planar_imaging, __version__, winston_lutz, utilities, watcher


class PylinacGUI(Frame):

    def __init__(self, master=None):
        super().__init__(master)
        self.pack()
        self.notebook = Notebook(self)
        self.init_pf()
        self.init_vmat()
        self.init_catphan()
        self.init_logs()
        # self.init_tg51()
        self.init_star()
        self.init_planar_imaging()
        self.init_winstonlutz()
        self.init_watcher()
        self.init_help()
        for child in self.winfo_children():
            child.grid_configure(padx=10, pady=10)

    def init_help(self):

        def upload():
            webbrowser.open('https://www.dropbox.com/request/YKRu4AmuPsXu55uQq761')

        def gotoforum():
            webbrowser.open('https://groups.google.com/forum/#!forum/pylinac')

        def gotogithub():
            webbrowser.open('https://github.com/jrkerns/pylinac')

        def gotortd():
            webbrowser.open('https://pylinac.readthedocs.io/en/latest/')

        self.help_tab = Frame(self.notebook)
        Label(self.help_tab, text='Having trouble?\nWant to donate your images for program improvement?\nUpload them below')
        Button(self.help_tab, text='Upload Files', command=upload)
        Label(self.help_tab, text='Complete documentation is available on ReadTheDocs')
        Button(self.help_tab, text='ReadTheDocs', command=gotortd)
        Label(self.help_tab, text='Help is also available in the Pylinac forum')
        Button(self.help_tab, text='Go to forum', command=gotoforum)
        Label(self.help_tab, text='The source code of this program and all analyses is available on Github')
        Button(self.help_tab, text='Github', command=gotogithub)
        self.notebook.add(self.help_tab, text='Help/Upload Images')
        for child in self.help_tab.winfo_children():
            child.grid_configure(padx=10, pady=5)

    def init_watcher(self):

        def load_yaml():
            f = filedialog.askopenfilename()
            self.watch_yaml.set(f)

        def load_dir():
            f = filedialog.askdirectory()
            self.watch_dir.set(f)

        def process_dir():
            watcher.process(directory=self.watch_dir.get(), config_file=self.watch_yaml.get(), force=self.watch_force.get())

        def open_dir():
            webbrowser.open(self.watch_dir.get())

        self.watch_tab = Frame(self.notebook)
        self.watch_yaml = StringVar(value=osp.join(osp.dirname(__file__), 'watcher_config.yml'))
        self.watch_dir = StringVar()
        self.watch_force = BooleanVar(value=False)
        Button(self.watch_tab, text='Load YAML config file...', command=load_yaml).grid(column=1, row=1)
        Label(self.watch_tab, textvariable=self.watch_yaml).grid(column=1, row=2)
        Button(self.watch_tab, text='Select analysis directory...', command=load_dir).grid(column=1, row=3)
        Label(self.watch_tab, textvariable=self.watch_dir).grid(column=1, row=4)
        Checkbutton(self.watch_tab, text='Force analysis (even previously analysed files)?', variable=self.watch_force).grid(column=1, row=5)
        Button(self.watch_tab, text='Process directory', command=process_dir).grid(column=1, row=6)
        Button(self.watch_tab, text='Open analysis directory', command=open_dir).grid(column=1, row=7)
        self.notebook.add(self.watch_tab, text='Batch Processor')
        for child in self.watch_tab.winfo_children():
            child.grid_configure(padx=10, pady=5)

    def init_vmat(self):

        def load_open():
            f = filedialog.askopenfilename()
            self.vmat_openimg.set(f)

        def load_dmlc():
            f = filedialog.askopenfilename()
            self.vmat_dmlcimg.set(f)

        def analyze_vmat():
            images = (self.vmat_openimg.get(), self.vmat_dmlcimg.get())
            if self.vmat_test.get() == 'DRGS':
                v = vmat.DRGS(image_paths=images)
            else:
                v = vmat.DRMLC(image_paths=images)
            v.analyze(tolerance=self.vmat_tol.get())
            fname = osp.join(self.vmat_dmlcimg.get().replace('.dcm', '.pdf'))
            fname = utilities.file_exists(fname)
            v.publish_pdf(fname)
            self.vmat_pdf.set(fname)
            utilities.open_path(fname)

        self.vmat_tab = Frame(self.notebook)
        self.vmat_openimg = StringVar()
        self.vmat_dmlcimg = StringVar()
        self.vmat_test = StringVar(value='DRGS')
        self.vmat_tol = DoubleVar(value=1.5)
        self.vmat_pdf = StringVar()
        Button(self.vmat_tab, text='Load Open Image...', command=load_open).grid(column=1, row=1)
        Button(self.vmat_tab, text='Load DMLC Image...', command=load_dmlc).grid(column=1, row=3)
        Label(self.vmat_tab, textvariable=self.vmat_openimg).grid(column=1, row=2)
        Label(self.vmat_tab, textvariable=self.vmat_dmlcimg).grid(column=1, row=4)
        Label(self.vmat_tab, text='Test type:').grid(column=1, row=5)
        Combobox(self.vmat_tab, values=('DRGS', 'DRMLC'), textvariable=self.vmat_test).grid(column=2, row=5)
        Label(self.vmat_tab, text='Tolerance (%):').grid(column=1, row=6)
        Entry(self.vmat_tab, width=7, textvariable=self.vmat_tol).grid(column=2, row=6)
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
                pf = picketfence.PicketFence(self.pf_file.get(), filter=3)
            else:
                pf = picketfence.PicketFence(self.pf_file.get())
            atol = self.pf_atol.get() if self.pf_atol.get() == 0 else None
            pickets = self.pf_pickets.get() if self.pf_pickets.get() == 0 else None
            hd = self.pf_hdmlc.get()
            pf.analyze(tolerance=self.pf_tol.get(),
                       action_tolerance=atol,
                       hdmlc=hd,
                       num_pickets=pickets,
                       )
            fname = osp.join(self.pf_file.get().replace('.dcm', '.pdf'))
            fname = utilities.file_exists(fname)
            pf.publish_pdf(fname)
            self.pf_pdf.set(fname)
            utilities.open_path(fname)

        self.pf_tab = Frame(self.notebook)
        self.pf_filter = BooleanVar(value=False)
        self.pf_file = StringVar()
        self.pf_tol = DoubleVar(value=0.5)
        self.pf_atol = DoubleVar(value=0.25)
        self.pf_pickets = IntVar(value=10)
        self.pf_hdmlc = BooleanVar(value=False)
        self.pf_pdf = StringVar()
        Checkbutton(self.pf_tab, text='Apply median filter', variable=self.pf_filter).grid(column=1, row=3)
        Button(self.pf_tab, text='Load File...', command=load_file).grid(column=1, row=1)
        Label(self.pf_tab, text='File:').grid(column=1, row=2)
        Label(self.pf_tab, textvariable=self.pf_file).grid(column=2, row=3)
        Label(self.pf_tab, text='Tolerance (mm):').grid(column=1, row=4)
        Entry(self.pf_tab, width=7, textvariable=self.pf_tol).grid(column=2, row=4)
        Label(self.pf_tab, text='Action Tolerance (mm):').grid(column=1, row=5)
        Entry(self.pf_tab, width=7, textvariable=self.pf_atol).grid(column=2, row=5)
        Label(self.pf_tab, text='Number of pickets:').grid(column=1, row=6)
        Entry(self.pf_tab, width=7, textvariable=self.pf_pickets).grid(column=2, row=6)
        Checkbutton(self.pf_tab, text='HD-MLC?', variable=self.pf_hdmlc).grid(column=1, row=7)
        Button(self.pf_tab, text='Analyze', command=analyze_pf).grid(column=1, row=8)
        Label(self.pf_tab, text='Analysis will analyze the file according to the settings, \nsave a PDF in the same directory as the original file location and then open it.').grid(column=1, row=9)
        self.notebook.add(self.pf_tab, text='Picket Fence')
        for child in self.pf_tab.winfo_children():
            child.grid_configure(padx=10, pady=5)

    def init_catphan(self):

        def load_dir():
            f = filedialog.askdirectory()
            self.ct_file.set(f)

        def load_zip():
            f = filedialog.askopenfilename()
            self.ct_file.set(f)

        def analyze_cbct():
            if osp.isdir(self.ct_file.get()):
                cat = getattr(ct, self.ct_catphantype.get())(self.ct_file.get())
                fname = osp.join(self.ct_file.get(), 'CBCT Analysis.pdf')
            else:
                cat = getattr(ct, self.ct_catphantype.get()).from_zip(self.ct_file.get())
                fname = self.ct_file.get().replace('.zip', '.pdf')
            cat.analyze(hu_tolerance=self.ct_hu.get(), thickness_tolerance=self.ct_thickness.get(),
                        scaling_tolerance=self.ct_scaling.get())
            fname = utilities.file_exists(fname)
            cat.publish_pdf(fname)
            self.ct_pdf.set(fname)
            utilities.open_path(fname)

        self.ct_tab = Frame(self.notebook)
        self.ct_file = StringVar()
        self.ct_catphantype = StringVar()
        self.ct_hu = IntVar(value=40)
        self.ct_scaling = DoubleVar(value=1)
        self.ct_thickness = DoubleVar(value=0.2)
        self.ct_pdf = StringVar()
        Label(self.ct_tab, text='Load EITHER a directory or ZIP file').grid(column=2, row=1)
        Button(self.ct_tab, text='Load Directory...', command=load_dir).grid(column=2, row=2)
        Button(self.ct_tab, text='Load ZIP file...', command=load_zip).grid(column=2, row=3)
        Label(self.ct_tab, textvariable=self.ct_file).grid(column=2, row=4)
        Label(self.ct_tab, text='CatPhan type:').grid(column=2, row=5)
        Combobox(self.ct_tab, values=('CatPhan504', 'CatPhan503', 'CatPhan600', 'CatPhan604'), textvariable=self.ct_catphantype).grid(column=2, row=6)
        Label(self.ct_tab, text='HU Tolerance (HU):').grid(column=1, row=7)
        Entry(self.ct_tab, width=7, textvariable=self.ct_hu).grid(column=1, row=8)
        Label(self.ct_tab, text='Scaling tolerance (mm):').grid(column=2, row=7)
        Entry(self.ct_tab, width=7, textvariable=self.ct_scaling).grid(column=2, row=8)
        Label(self.ct_tab, text='Thickness tolerance (mm):').grid(column=3, row=7)
        Entry(self.ct_tab, width=7, textvariable=self.ct_thickness).grid(column=3, row=8)
        Button(self.ct_tab, text='Analyze', command=analyze_cbct).grid(column=2, row=9)
        Label(self.ct_tab,
              text='Analysis will analyze the file(s) according to the settings, \nsave a PDF in the same directory as the original file location and then open it.').grid(
            column=2, row=10)
        Label(self.ct_tab, text='Save file:').grid(column=2, row=11)
        Label(self.ct_tab, textvariable=self.ct_pdf).grid(column=2, row=12)
        self.notebook.add(self.ct_tab, text='CatPhan')
        for child in self.ct_tab.winfo_children():
            child.grid_configure(padx=10, pady=5)

    def init_logs(self):

        def load_log():
            f = filedialog.askopenfilename()
            self.log_file.set(f)

        def analyze_log():
            log = log_analyzer.load_log(self.log_file.get())
            name, _ = osp.splitext(self.log_file.get())
            fname = name + '.pdf'
            fname = utilities.file_exists(fname)
            log.publish_pdf(fname)
            self.log_pdf.set(fname)
            utilities.open_path(fname)

        self.log_tab = Frame(self.notebook)
        self.log_file = StringVar()
        self.log_pdf = StringVar()
        Button(self.log_tab, text='Load log file...', command=load_log).grid(column=1, row=1)
        Label(self.log_tab, textvariable=self.log_file).grid(column=1, row=2)
        Button(self.log_tab, text='Analyze', command=analyze_log).grid(column=1, row=9)
        Label(self.log_tab,
              text='Analysis will analyze the file(s) according to the settings, \nsave a PDF in the same directory as the original file location and then open it.').grid(
            column=1, row=10)
        Label(self.log_tab, text='Save file:').grid(column=1, row=11)
        Label(self.log_tab, textvariable=self.log_pdf).grid(column=1, row=12)
        self.notebook.add(self.log_tab, text='Machine Logs')
        for child in self.log_tab.winfo_children():
            child.grid_configure(padx=10, pady=5)

    def init_star(self):

        def load_star():
            f = filedialog.askopenfilename()
            self.star_file.set(f)

        def analyze_star():
            star = starshot.Starshot(self.star_file.get(), sid=self.star_sid.get(),
                                     dpi=self.star_dpi.get())
            star.analyze(radius=self.star_radius.get(), tolerance=self.star_tolerance.get(),
                         recursive=self.star_recursive.get())
            name, _ = osp.splitext(self.star_file.get())
            fname = name + '.pdf'
            fname = utilities.file_exists(fname)
            star.publish_pdf(fname)
            self.star_pdf.set(fname)
            utilities.open_path(fname)

        self.star_tab = Frame(self.notebook)
        self.star_file = StringVar()
        self.star_pdf = StringVar()
        self.star_dpi = DoubleVar()
        self.star_sid = DoubleVar()
        self.star_radius = DoubleVar(value=0.85)
        self.star_tolerance = DoubleVar(value=1)
        self.star_recursive = BooleanVar(value=True)
        Button(self.star_tab, text='Load starshot file...', command=load_star).grid(column=1, row=1)
        Label(self.star_tab, textvariable=self.star_file).grid(column=1, row=2)
        Label(self.star_tab, text='DPI (if file is not DICOM):').grid(column=1, row=3)
        Entry(self.star_tab, width=7, textvariable=self.star_dpi).grid(column=2, row=3)
        Label(self.star_tab, text='SID (mm; if file is not DICOM):').grid(column=1, row=4)
        Entry(self.star_tab, width=7, textvariable=self.star_sid).grid(column=2, row=4)
        Label(self.star_tab, text='Normalized analysis radius (0.2-1.0):').grid(column=1, row=5)
        Entry(self.star_tab, width=7, textvariable=self.star_radius).grid(column=2, row=5)
        Checkbutton(self.star_tab, text='Recursive analysis?', variable=self.star_recursive).grid(column=1, row=6)
        Label(self.star_tab, text='Tolerance (mm):').grid(column=1, row=7)
        Entry(self.star_tab, width=7, textvariable=self.star_tolerance).grid(column=2, row=7)
        Button(self.star_tab, text='Analyze', command=analyze_star).grid(column=1, row=9)
        Label(self.star_tab,
              text='Analysis will analyze the file(s) according to the settings, \nsave a PDF in the same directory as the original file location and then open it.').grid(
              column=1, row=10)
        Label(self.star_tab, text='Save file:').grid(column=1, row=11)
        Label(self.star_tab, textvariable=self.star_pdf).grid(column=1, row=12)
        self.notebook.add(self.star_tab, text='Starshot')
        for child in self.star_tab.winfo_children():
            child.grid_configure(padx=10, pady=5)

    def init_planar_imaging(self):

        def load_phan():
            f = filedialog.askopenfilename()
            self.phan_file.set(f)

        def analyze_phan():
            phantom = getattr(planar_imaging, self.phan_type.get())(self.phan_file.get())
            phantom.analyze()
            name, _ = osp.splitext(self.phan_file.get())
            fname = name + '.pdf'
            fname = utilities.file_exists(fname)
            phantom.publish_pdf(utilities.file_exists(fname))
            self.phan_pdf.set(fname)
            utilities.open_path(fname)

        self.phan_tab = Frame(self.notebook)
        self.phan_file = StringVar()
        self.phan_pdf = StringVar()
        self.phan_locon = DoubleVar(value=0.1)
        self.phan_hicon = DoubleVar(value=0.5)
        self.phan_inver = BooleanVar(value=False)
        self.phan_type = StringVar(value='LeedsTOR')
        Button(self.phan_tab, text='Load planar phantom DICOM file...', command=load_phan).grid(column=1, row=1)
        Label(self.phan_tab, textvariable=self.phan_file).grid(column=1, row=2)
        Label(self.phan_tab, text='Phantom:').grid(column=1, row=3)
        Combobox(self.phan_tab, values=('LeedsTOR', 'LasVegas', 'StandardImagingQC3'), textvariable=self.phan_type).grid(column=2, row=3)
        Label(self.phan_tab, text='Low contrast threshold:').grid(column=1, row=4)
        Entry(self.phan_tab, width=7, textvariable=self.phan_locon).grid(column=2, row=4)
        Label(self.phan_tab, text='High contrast threshold:').grid(column=1, row=5)
        Entry(self.phan_tab, width=7, textvariable=self.phan_hicon).grid(column=2, row=5)
        Checkbutton(self.phan_tab, text='Force image inversion?', variable=self.phan_inver).grid(column=1, row=6)
        Button(self.phan_tab, text='Analyze', command=analyze_phan).grid(column=1, row=9)
        Label(self.phan_tab,
              text='Analysis will analyze the file(s) according to the settings, \nsave a PDF in the same directory as the original file location and then open it.').grid(
              column=1, row=10)
        Label(self.phan_tab, text='Save file:').grid(column=1, row=11)
        Label(self.phan_tab, textvariable=self.phan_pdf).grid(column=1, row=12)
        self.notebook.add(self.phan_tab, text='2D Phantoms')
        for child in self.phan_tab.winfo_children():
            child.grid_configure(padx=10, pady=5)

    def init_winstonlutz(self):

        def load_dir():
            f = filedialog.askdirectory()
            self.wl_file.set(f)

        def load_zip():
            f = filedialog.askopenfilename()
            self.wl_file.set(f)

        def analyze_wl():
            if osp.isdir(self.wl_file.get()):
                wl = winston_lutz.WinstonLutz(self.wl_file.get())
                fname = osp.join(self.wl_file.get(), 'W-L Analysis.pdf')
            else:
                wl = winston_lutz.WinstonLutz.from_zip(self.wl_file.get())
                fname = self.wl_file.get().replace('.zip', '.pdf')
            fname = utilities.file_exists(fname)
            wl.publish_pdf(fname)
            self.wl_pdf.set(fname)
            utilities.open_path(fname)

        self.wl_tab = Frame(self.notebook)
        self.wl_file = StringVar()
        self.wl_pdf = StringVar()
        Label(self.wl_tab, text='Load EITHER a directory or ZIP file').grid(column=2, row=1)
        Button(self.wl_tab, text='Load Directory...', command=load_dir).grid(column=2, row=2)
        Button(self.wl_tab, text='Load ZIP file...', command=load_zip).grid(column=2, row=3)
        Label(self.wl_tab, textvariable=self.wl_file).grid(column=2, row=4)
        Button(self.wl_tab, text='Analyze', command=analyze_wl).grid(column=2, row=9)
        Label(self.wl_tab,
              text='Analysis will analyze the file(s) according to the settings, \nsave a PDF in the same directory as the original file location and then open it.').grid(
            column=2, row=10)
        Label(self.wl_tab, text='Save file:').grid(column=2, row=11)
        Label(self.wl_tab, textvariable=self.wl_pdf).grid(column=2, row=12)
        self.notebook.add(self.wl_tab, text='Winston-Lutz')
        for child in self.wl_tab.winfo_children():
            child.grid_configure(padx=10, pady=5)

    def init_tg51(self):

        self.tg_tab = Frame(self.notebook)
        self.tg_pdf = StringVar()
        self.tg_temp = DoubleVar(value=22)
        self.tg_press = DoubleVar(value=760)
        r, r2 = itertools.count(), itertools.count()
        Label(self.tg_tab, text='Temperature (C):').grid(column=1, row=next(r))
        Entry(self.tg_tab, width=7, textvariable=self.tg_temp).grid(column=2, row=next(r2))
        Label(self.tg_tab, text='Pressure (mmHg):').grid(column=1, row=next(r))
        Entry(self.tg_tab, width=7, textvariable=self.tg_press).grid(column=2, row=next(r2))
        self.notebook.add(self.tg_tab, text='TG-51')
        for child in self.tg_tab.winfo_children():
            child.grid_configure(padx=10, pady=5)


def gui():

    def on_exit():
        if messagebox.askokcancel("Quit", "Do you want to quit?"):
            root.quit()

    root = Tk()
    root.title('Pylinac GUI ' + __version__)
    root.protocol("WM_DELETE_WINDOW", on_exit)
    app = PylinacGUI(master=root)
    app.mainloop()
    root.destroy()
    del root
