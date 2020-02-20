# import standard or third party modules
import os
import re
import tkinter as tk
import traceback
from shutil import which
from subprocess import Popen
from time import time
from tkinter import ttk, filedialog, messagebox, font

from source.configuration import cfg, gui, adjust_dir_path, adjust_file_path, find_exe
from source.converter import MutationVcfMerger, UniProtConverter, PatricConverter, RegulonDBConverter
from source.database import Database
from source.log import ALog
# import own modules
from source.stats import Stats
from source.tools import check_dir, setup_directory, clean_directory

from source.ngs import NGSPipeline

__author__ = 'Markus Hollander'
__copyright__ = 'Copyright (C) 2017 Markus Hollander'
__license__ = 'GPL V3'
__version__ = '2.0'


def out_warn(out_dir):
    """
    Warn the user that the old output is going to be deleted.
    """
    return messagebox.askokcancel('Warning',
                                  'This is going to delete all files in the specified output directory.\n\n'
                                  'Directory path: {0}\n\n'
                                  'Press OK to continue, or Cancel to abort.'.format(out_dir))


class Key:
    """ Can be used to give other classes support for mouse-over. """
    def __init__(self, key):
        self.key = key
        self.label = gui.desc[key].label
        self.desc = gui.desc[key].desc
        
    def mouse_over(self, desc):
        """
        Activates mouse-over behaviour for the widget.
        :param desc: description box
        """
        self.bind('<Enter>', lambda _: desc.show(self.label, self.desc))
        self.bind('<Leave>', desc.clear)


class TLabel(ttk.Label, Key):
    """ A title label for sections_verbose that supports mouse-over. """
    def __init__(self, master, key):
        Key.__init__(self, key)
        ttk.Label.__init__(self, master, text=self.label, style=gui.st_style)

    def place(self, row):
        self.grid(row=row, column=0, columnspan=5, pady=gui.st_pady, padx=gui.st_padx, sticky=tk.W)


class STLabel(ttk.Label, Key):
    """ A title label for subsections that supports mouse-over. """
    def __init__(self, master, key):
        Key.__init__(self, key)
        ttk.Label.__init__(self, master, text=self.label, style=gui.sst_style)

    def place(self, row):
        self.grid(row=row, column=1, pady=gui.sst_pady, padx=gui.sst_padx, sticky=tk.W)


class SLabel(ttk.Label, Key):
    """ A status label that supports mouse-over. """
    def __init__(self, master, key, var, pos='l'):
        Key.__init__(self, key)
        ttk.Label.__init__(self, master, textvariable=var, style=gui.sl_style)
        if pos == 'c':
            self.configure(anchor=tk.CENTER, width=gui.b_width)
        self.pos = pos

    def place(self, row):
        if self.pos == 'l':
            self.grid(row=row, column=0, sticky=tk.W)
        elif self.pos == 'r':
            self.grid(row=row, column=1, sticky=tk.E)
        else:
            self.grid(row=row, column=0, columnspan=2, sticky=tk.EW)


class ELabel(ttk.Label, Key):
    """ A label naming an entry field that supports mouse-over. """
    def __init__(self, master, key):
        Key.__init__(self, key)
        ttk.Label.__init__(self, master, text=self.label + ':', style=gui.el_style)

    def place(self, row):
        self.grid(row=row, column=1, pady=gui.el_pady, padx=gui.el_padx, sticky=tk.W)


class RButton(ttk.Button, Key):
    """ Main run button that supports mouse-over. """
    def __init__(self, master, key, cmd):
        Key.__init__(self, key)
        ttk.Button.__init__(self, master, text=self.label, command=cmd, width=gui.rb_width, style=gui.rb_style)

    def place(self, row):
        self.grid(row=row, column=3, columnspan=2, sticky=tk.E, pady=gui.rb_pady, padx=gui.rb_padx)


class LSButton(ttk.Button, Key):
    """ Setting window left button that supports mouse-over. """
    def __init__(self, master, key, cmd):
        Key.__init__(self, key)
        ttk.Button.__init__(self, master, text=self.label, command=cmd, width=gui.lsb_width, style=gui.rb_style)

    def place(self, row):
        self.grid(row=row, column=2, sticky=tk.E, pady=gui.lsb_pady, padx=gui.lsb_padx)


class RSButton(ttk.Button, Key):
    """ Setting window right button that supports mouse-over. """
    def __init__(self, master, key, cmd):
        Key.__init__(self, key)
        ttk.Button.__init__(self, master, text=self.label, command=cmd, width=gui.rsb_width, style=gui.rb_style)

    def place(self, row):
        self.grid(row=row, column=3, columnspan=2, sticky=tk.E, pady=gui.rsb_pady, padx=gui.rsb_padx)


class FButton(ttk.Button, Key):
    """ File/directory chooser button that supports mouse-over. """
    def __init__(self, master, key, cmd):
        Key.__init__(self, key)
        ttk.Button.__init__(self, master, text='...', command=cmd, style=gui.fb_style, width=gui.fb_width)


class CEntry(ttk.Entry, Key):
    """ Custom entry field that supports mouse-over. """
    def __init__(self, master, key, var, width):
        Key.__init__(self, key)
        ttk.Entry.__init__(self, master, textvariable=var, width=width, font=gui.ef_font)


class TCheckbutton(ttk.Checkbutton, Key):
    """ Checkbutton with label that supports mouse-over. """
    def __init__(self, master, key, default):
        Key.__init__(self, key)
        self.var = tk.BooleanVar(value=default)
        ttk.Checkbutton.__init__(self, master, text=self.label, variable=self.var, style=gui.cb_style)

    def place(self, row):
        self.grid(row=row, column=1, columnspan=3, stick=tk.W, pady=gui.tcb_pady, padx=gui.tcb_padx)


class ECheckbutton(ttk.Checkbutton, Key):
    """ Checkbutton without label that supports mouse-over and takes a command. """
    def __init__(self, master, key, var, cmd):
        Key.__init__(self, key)
        if cmd:
            ttk.Checkbutton.__init__(self, master, text='', variable=var, command=cmd)
        else:
            ttk.Checkbutton.__init__(self, master, text='', variable=var)

    def place(self, row):
        self.grid(row=row, column=0, sticky=tk.W, pady=gui.ecb_pady, padx=gui.ecb_padx)


class Sep(ttk.Separator, Key):
    def __init__(self, master, key):
        Key.__init__(self, key)
        ttk.Separator.__init__(self, master, orient=tk.HORIZONTAL)


""" AGGREGATES """


class TSep(ttk.Frame):
    """ Frame that contains an (optional) title and a separator line that supports mouse-over. """
    def __init__(self, master, key):
        ttk.Frame.__init__(self, master)
        self.columnconfigure(1, weight=1)

        if gui.desc[key]:
            self.title = STLabel(self, key)
        else:
            self.title = None

        self.sep = Sep(self, key)

    def place(self, row):
        if self.title:
            self.title.grid(row=0, column=0, sticky=tk.W, pady=gui.sst_pady, padx=gui.sst_padx)
        self.sep.grid(row=0, column=1, sticky=tk.EW, pady=cfg.gui.tss_pady, padx=gui.tss_padx)
        self.grid(row=row, column=1, columnspan=4, sticky=tk.EW)

    def mouse_over(self, desc):
        """
        Activates mouse-over behaviour.
        """
        if self.title:
            self.title.mouse_over(desc)
        self.sep.mouse_over(desc)


class NInput:
    """
    Input line for numbers that supports mouse-over: title - entry field.
    """
    def __init__(self, master, key, default):
        self.key = key
        self.title = ELabel(master, key)
        self.var = tk.StringVar(value=default)
        self.field = CEntry(master, key, self.var, gui.ne_width)

    def show(self, row):
        """
        Shows the input fields in the application window.
        """
        self.title.grid(row=row, column=1, sticky=tk.W, pady=gui.el_pady, padx=gui.el_padx)
        self.field.grid(row=row, column=2, sticky=tk.W, pady=gui.ef_pady, padx=gui.ef_padx, columnspan=2)

    def hide(self):
        """
        Removes the input fields from the application window.
        """
        self.title.grid_remove()
        self.field.grid_remove()

    def mouse_over(self, desc):
        """
        Activates mouse-over behaviour.
        """
        self.title.mouse_over(desc)
        self.field.mouse_over(desc)

    def validate(self):
        pass


class PercInput(NInput):
    """
    Input line for percentages that supports mouse-over: title - entry field.
    """
    def __init__(self, master, key, default):
        NInput.__init__(self, master, key, default)

    def validate(self):
        """
        Checks if value is a valid percentage.
        :return: True if the number is a percentage, False otherwise
        """
        try:
            val = float(self.var.get())

            if val < 0 or val > 1:
                messagebox.showerror('Input Error', '{0} needs to be between 0 and 1.'.format(gui.desc[self.key].label))
                return False

            return True

        except ValueError:
            messagebox.showerror('Input Error', '{0} needs to be a number.'.format(gui.desc[self.key].label))
            return False


class IntInput(NInput):
    """
    Input line for integers that supports mouse-over: title - entry field.
    """
    def __init__(self, master, key, default, pos=False, neg=False):
        NInput.__init__(self, master, key, default)
        self.pos = pos
        self.neg = neg

    def validate(self):
        """
        Checks if value is a valid number.
        :return: True if the number is valid, False otherwise
        """
        try:
            val = int(self.var.get())

            if self.pos and val < 0:
                messagebox.showerror('Input Error', '{0} needs to be >= 0.'.format(gui.desc[self.key].label))
                return False

            if self.neg and val > 0:
                messagebox.showerror('Input Error', '{0} needs to be <= 0.'.format(gui.desc[self.key].label))
                return False

            return True

        except ValueError:
            messagebox.showerror('Input Error', '{0} needs to be a number.'.format(gui.desc[self.key].label))
            return False


class TextInput:
    """
    Input line for text.
    """
    def __init__(self, master, key, default):
        self.key = key
        self.title = ELabel(master, key)
        self.var = tk.StringVar(value=default)
        self.field = CEntry(master, key, self.var, gui.fe_width)

    def hide(self):
        """
        Removes the input field from the application window.
        """
        self.title.grid_remove()
        self.field.grid_remove()

    def show(self, row, label=True):
        """
        Shows the input fields in the application window.
        """
        if label:
            self.title.grid(row=row, column=1, sticky=tk.W, pady=gui.el_pady, padx=gui.el_padx)
        self.field.grid(row=row, column=2, columnspan=2, sticky=tk.W, pady=gui.ef_pady, padx=gui.ef_padx)

    def mouse_over(self, desc):
        """
        Activates mouse-over behaviour.
        """
        self.title.mouse_over(desc)
        self.field.mouse_over(desc)

    def validate(self):
        """
        Checks if the input only contains allowed characters.
        :return: True if input is valid, False otherwise
        """
        self.var.set(' '.join(self.var.get().split()).strip())

        if not re.match(r'[\w -]*$', self.var.get()):
            messagebox.showerror('Input Error', '{0} may only contain letters, numbers, underscores, dashes and '
                                                'spaces.'.format(gui.desc[self.key].label))
            return False
        return True


class FDInput:
    """
    Input line for file/directory paths that supports mouse-over: title - entry field - chooser button.
    """
    def __init__(self, master, key, default):
        self.title = ELabel(master, key)
        self.var = tk.StringVar(value=default)
        self.field = CEntry(master, key, self.var, gui.fe_width)
        self.chooser = FButton(master, key, self.ask)
        self.header = gui.desc[key].label
        self.titlefy()

    def hide(self):
        """
        Removes the input field from the application window.
        """
        self.title.grid_remove()
        self.field.grid_remove()
        self.chooser.grid_remove()

    def show(self, row, label=True):
        """
        Shows the input fields in the application window.
        """
        if label:
            self.title.grid(row=row, column=1, sticky=tk.W, pady=gui.el_pady, padx=gui.el_padx)
        self.field.grid(row=row, column=2, columnspan=2, sticky=tk.W, pady=gui.ef_pady, padx=gui.ef_padx)
        self.chooser.grid(row=row, column=4, sticky=tk.E, pady=gui.fb_pady, padx=gui.fb_padx)

    def mouse_over(self, desc):
        """
        Activates mouse-over behaviour.
        """
        self.title.mouse_over(desc)
        self.field.mouse_over(desc)
        self.chooser.mouse_over(desc)

    def titlefy(self):
        """
        Transforms the name to title case.
        """
        words = self.header.split()

        for i in range(len(words)):
            if not words[i].isupper():
                words[i] = words[i].title()

        self.header = ' '.join(words)

    def ask(self):
        pass

    def validate(self):
        pass


class DirInput(FDInput):
    """
    Input field for a directory.
    """
    def __init__(self, master, key, default):
        FDInput.__init__(self, master, key, default)
        self.name = gui.desc[key].label.lower().replace('directory', '').strip()

    def ask(self):
        """
        Open the directory chooser dialog.
        """
        p = filedialog.askdirectory(initialdir=(os.path.abspath(self.var.get())),
                                    title='Select {0}'.format(self.header))
        if p:
            self.var.set(adjust_dir_path(p))

    def open(self):
        """
        Opens the directory in the OS file explorer.
        """
        path = os.path.abspath(self.var.get())

        # open command depending on the operating system
        if cfg.os == cfg.linux:
            Popen(["xdg-open", path])
        elif cfg.os == cfg.mac:
            Popen(["open", path])
        elif cfg.os == cfg.windows:
            os.startfile(path)
        elif cfg.os == cfg.cygwin:
            Popen(["xdg-open", path])

    def validate(self):
        """
        Checks if the directory exists and is a directory, and shows an error message if it is not.
        :return: True if the directory is valid, False otherwise
        """
        return check_dir(self.name, self.var.get())

    def setup(self):
        """
        Sets up the specified directory if it does not exist yet and shows and error message if that fails.
        :return: True if the directory already existed/could be created and is accessible, False otherwise
        """
        return setup_directory(self.name, self.var.get())

    def clean(self):
        """
        Deletes all files and subdirectories unless they are exempt.
        :return: True if all files/subdirectories could be delete, False otherwise
        """
        return clean_directory(self.name, self.var.get(), True)


class ReadFileInput(FDInput):
    """
    Input field for files.
    """
    def __init__(self, master, key, default):
        FDInput.__init__(self, master, key, default)
        self.name = gui.desc[key].label.lower().replace('file', '').strip()

    def ask(self):
        """
        Open the file chooser dialog.
        """
        d, f = os.path.split(os.path.abspath(self.var.get()))
        p = filedialog.askopenfilename(initialfile=f,
                                       initialdir=d,
                                       title='Select {0}'.format(self.header))
        if p:
            self.var.set(adjust_file_path(p))

    def validate(self):
        """
        Checks if the file exists, is a file and if it is readable and shows an error message if it is not.
        :return: True if the file is valid, False otherwise
        """
        path = adjust_file_path(self.var.get())

        if not os.path.isfile(path):
            messagebox.showerror('File Error', 'The {0} file does not exist.\n\n'
                                               'File: {1}'.format(self.name, path))
            return False
        if not os.access(path, os.R_OK):
            messagebox.showerror('File Error', 'The {0} file is not readable.\n\n'
                                               'File {1}'.format(self.name, path))
            return False
        return True

    def open(self):
        """
        Open the directory containing the file.
        """
        path = os.path.dirname(os.path.abspath(self.var.get()))

        # open command depending on the operating system
        if cfg.os == cfg.linux:
            Popen(["xdg-open", path])
        elif cfg.os == cfg.mac:
            Popen(["open", path])
        elif cfg.os == cfg.windows:
            os.startfile(path)
        elif cfg.os == cfg.cygwin:
            Popen(["xdg-open", path])


class SaveFileInput(FDInput):
    """
    Input field for files.
    """
    def __init__(self, master, key, default):
        FDInput.__init__(self, master, key, default)
        self.name = gui.desc[key].label.lower().replace('file', '').strip()

    def ask(self):
        """
        Open the file chooser dialog.
        """
        d, f = os.path.split(os.path.abspath(self.var.get()))
        p = filedialog.asksaveasfilename(initialfile=f,
                                         initialdir=d,
                                         title='Select {0}'.format(self.header))
        if p:
            self.var.set(adjust_file_path(p))

    def validate(self, exists=True):
        """
        Checks if the file exists, is a file and if it is readable and shows an error message if it is not.
        :return: True if the file is valid, False otherwise
        """
        path = adjust_file_path(self.var.get())
        dir = '/'.join(path.split('/')[:-1])

        if not os.path.exists(dir):
            os.makedirs(dir)

        if exists:
            if not os.path.isfile(path):
                messagebox.showerror('File Error', 'The {0} file does not exist.\n\n'
                                                   'File: {1}'.format(self.name, path))
                return False
            if not os.access(path, os.W_OK):
                messagebox.showerror('File Error', 'The {0} file is not writable.\n\n'
                                                   'File {1}'.format(self.name, path))
                return False
        if not os.path.isfile(path):
            open(path, 'w').close()
        return True

    def open(self):
        """
        Open the directory containing the file.
        """
        path = os.path.dirname(os.path.abspath(self.var.get()))

        # open command depending on the operating system
        if cfg.os == cfg.linux:
            Popen(["xdg-open", path])
        elif cfg.os == cfg.mac:
            Popen(["open", path])
        elif cfg.os == cfg.windows:
            os.startfile(path)
        elif cfg.os == cfg.cygwin:
            Popen(["xdg-open", path])


class ExeInput(FDInput):
    """
    Input field for executables.
    """
    def __init__(self, master, key, default, name):
        FDInput.__init__(self, master, key, default)
        self.name = name

    def ask(self):
        """
        Open the file chooser dialog.
        """
        self.var.set(adjust_file_path(filedialog.askopenfilename()))

    def validate(self):
        """
        Checks if a path points to a valid executable and shows an error message if not.
        :return: True if it is valid, False otherwise
        """
        path = find_exe(self.name, adjust_file_path(self.var.get()))

        if not which(path):
            msg = 'The {0} path does not point to an executable. Please select a valid executable ' \
                  'under \'Settings\' --> \'Advanced NGS settings\'.\n\nExecutable path: {1}. '
            messagebox.showerror('Executable Path Error',
                                 msg.format(self.name, path))
            return False

        self.var.set(adjust_file_path(path))

        return True


class Selector:
    """
    Check button with title and separator line, as well as associated input fields that can be toggled
    via the check button.
    """
    def __init__(self, master, key, default, inputs, label=True):
        self.row = 0
        self.var = tk.BooleanVar(value=default)
        self.check = ECheckbutton(master, key, self.var, None)
        self.sep = TSep(master, key)
        self.inputs = inputs
        self.label = label

    def mouse_over(self, desc):
        """
        Activates mouse-over behaviour.
        """
        self.check.mouse_over(desc)
        self.sep.mouse_over(desc)

        for input in self.inputs:
            input.mouse_over(desc)

    def place(self, row):
        """
        Places the selector and the associated input fields in the application window.
        """
        self.row = row
        self.check.place(row)
        self.sep.place(row)

        for i in range(len(self.inputs)):
            self.inputs[i].show(self.row + i + 1, self.label)


class Description(ttk.Frame):
    """
    A section showing descriptions of widgets.
    """
    def __init__(self, master, row):
        ttk.Frame.__init__(self, master, width=gui.dd_width)
        # default values
        self.title_text = gui.desc[gui.key.g_desc].label
        self.desc_text = gui.desc[gui.key.g_desc].desc
        # text variables
        self.title_var = tk.StringVar(value=self.title_text)
        self.desc_var = tk.StringVar(value=self.desc_text)
        # title and description
        self.title = ttk.Label(self, textvariable=self.title_var, style=gui.dt_style)
        self.desc = ttk.Label(self, textvariable=self.desc_var)
        self.desc.configure(font=gui.dd_font, wraplength=gui.dd_width, foreground=gui.dd_colour, anchor=tk.W)

        # place in the application window
        self.place(row)

    def place(self, row):
        """
        Places the section and its widgets in the window.
        """
        self.title.grid(row=0, column=0, sticky=tk.W, pady=gui.dt_pady, padx=gui.dt_padx)
        self.desc.grid(row=1, column=0, sticky=tk.EW, pady=gui.dd_pady, padx=gui.dd_padx)
        self.grid(row=row, column=0, columnspan=5, sticky=(tk.N, tk.W, tk.E))

    def show(self, title, desc):
        """
        Shows the specified title and description.
        """
        self.desc_var.set(desc)
        self.title_var.set(title)

    def clear(self, _):
        """
        Shows the default title and description.
        """
        self.desc_var.set(self.desc_text)
        self.title_var.set(self.title_text)


class Status:
    """
    A section in the main application window with status information of the currently running process.
    """
    def __init__(self, master, app, desc, row):
        """
        Initialises the frames and labels with their variables, places them in the application window and
        activates mouse-over.
        """
        self.rows = 11

        """ FRAMES """
        self.app = app
        outer = ttk.Frame(master)
        box = ttk.Frame(outer, relief=tk.SUNKEN, borderwidth=1)
        box.columnconfigure(0, weight=1)

        outer.grid(row=row, column=0, sticky=(tk.S, tk.E, tk.W))
        box.grid(row=1, column=0, columnspan=5, pady=gui.sf_pady, padx=gui.sf_padx, sticky=tk.EW)

        """ TEXT-VARIABLES """
        self.header = tk.StringVar()
        self.l_last = tk.StringVar()
        self.r_last = tk.StringVar()
        self.l_vars = dict()
        self.r_vars = dict()
        self.ind = '    '

        for i in range(self.rows):
            self.l_vars[i] = tk.StringVar()
            self.r_vars[i] = tk.StringVar()

        """ LABELS """
        key = gui.key.s_title
        # section title
        title = TLabel(outer, key)
        title.mouse_over(desc)
        title.place(row=0)
        # status box header
        header = SLabel(box, key, self.header, 'c')
        header.mouse_over(desc)
        header.place(0)
        # total running time
        l_last = SLabel(box, key, self.l_last, 'l')
        l_last.mouse_over(desc)
        l_last.place(14)
        r_last = SLabel(box, key, self.r_last, 'r')
        r_last.mouse_over(desc)
        r_last.place(14)
        # status lines
        for i in range(self.rows):
            left = SLabel(box, key, self.l_vars[i], 'l')
            left.mouse_over(desc)
            left.place(i + 2)
            right = SLabel(box, key, self.r_vars[i], 'r')
            right.mouse_over(desc)
            right.place(i + 2)

        """ SEPARATORS """
        Sep(box, key).grid(row=1, column=0, columnspan=2, sticky=tk.EW)
        Sep(box, key).grid(row=13, column=0, columnspan=2, sticky=tk.EW)

    def name(self, title):
        """
        Writes the name of the current process into the first line.
        """
        self.header.set(title)
        self.app.update()

    def running(self, line, text, indent=False):
        """
        Writes a in progress message into the specified line with the specified and optionally indented description.
        """
        text = self.text(text, indent)
            
        self.l_vars[line].set(text)
        self.r_vars[line].set('running')
        self.app.update()

    def fail(self, line, text, indent=False):
        """
        Writes a fail message into the specified line with the specified and optionally indented description.
        """
        text = self.text(text, indent)
        
        self.l_vars[line].set(text)
        self.r_vars[line].set('failed')
        self.app.update()

    def skip(self, line, text, indent=False):
        """
        Writes a skipp message into the specified line with the specified and optionally indented description.
        """
        text = self.text(text, indent)
        
        self.l_vars[line].set(text)
        self.r_vars[line].set('skipped')
        self.app.update()

    def time(self, line, text, t, indent=False):
        """
        Writes the running time line to the specified line with the specified and optionally indented description.
        """
        text = self.text(text, indent)
        
        self.l_vars[line].set(text)
        self.r_vars[line].set('{:.1f}s'.format(t))
        self.app.update()

    def of(self, line, text, left, right, indent=False):
        """
        Writes the running time line to the specified line with the specified and optionally indented description.
        """
        text = self.text(text, indent)

        self.l_vars[line].set(text)
        self.r_vars[line].set('{0} of {1}'.format(left, right))
        self.app.update()

    def total_time(self, t):
        """
        Write the total running time into the last line.
        :param t: time
        """
        self.l_last.set('total running time:')
        self.r_last.set('{:.1f}s'.format(t))
        self.app.update()

    def clear(self, start=0):
        """
        Remove previous status messages.
        """
        if start == 0:
            self.header.set('')
        self.l_last.set('')
        self.r_last.set('')

        for i in range(start, self.rows):
            self.r_vars[i].set('')
            self.l_vars[i].set('')

        self.app.update()

    def text(self, txt, indent=False):
        """
        Formats the name of the subroutine.
        """
        if indent:
            txt = self.ind + txt
        return txt + ':'


class NWindow(tk.Toplevel):
    """
    Window for additional settings or functionality.
    """
    def __init__(self, master):
        # build a new window that stays on top of the application window
        tk.Toplevel.__init__(self, master)
        self.title(gui.title)
        self.transient()
        self.grab_set()
        # set up the left and right content frames, as well as window borders
        self.left, self.right = self.setup()
        # set up the description section
        self.desc = Description(self.right, 0)
        self.desc.configure(height=gui.df_height, width=gui.df_width)
        self.desc.grid_propagate(False)
        # centre the window
        self.centre()
        
    def setup(self):
        """
        Sets up window borders and section frames.
        :return: left frame, right frame
        """
        # set up the window border frames and the separating frame
        ttk.Frame(self).grid(row=0, column=0, sticky=tk.NSEW, ipadx=gui.wl_padx)
        ttk.Frame(self).grid(row=0, column=4, sticky=tk.NSEW, ipadx=gui.wr_padx)
        ttk.Frame(self).grid(row=0, column=2, sticky=tk.NSEW, ipadx=gui.col_padx)
        # don't show bottom padding on Mac OS
        if cfg.os != cfg.mac:
            ttk.Frame(self).grid(row=1, column=0, ipady=gui.wb_pady)
        self.rowconfigure(0, weight=1)
        # set up the left and right section frame
        left = ttk.Frame(self)
        right = ttk.Frame(self)
        left.grid(row=0, column=1, sticky=tk.NSEW)
        right.grid(row=0, column=3, sticky=tk.NSEW)
        # change the expanding behaviour of the right frame
        right.rowconfigure(1, weight=1)
        left.rowconfigure(100, weight=1)

        return left, right

    def centre(self):
        """
        Places the application window in the middle of the screen.
        """
        self.update()
        # compute the x and y coordinates of the window using the screen width and height
        x = (self.winfo_screenwidth() - self.winfo_reqwidth()) // 2
        y = (self.winfo_screenheight() - self.winfo_reqheight()) // 2

        # set the window position (x coordinate, y coordinate) on the screen
        self.geometry('+{0}+{1}'.format(x, y))
        self.update()

    def cancel(self):
        """
        Discards changes and closes the window.
        """
        self.destroy()

    def seperator(self, row):
        """
        Places a separatur line in the specified row.
        """
        sep = ttk.Separator(self.left, orient=tk.HORIZONTAL)
        sep.grid(row=row, column=1, columnspan=4, pady=gui.sp_pady, padx=gui.sp_padx, sticky=tk.EW)

    def buttons(self, key, cmd, row):
        """
        Places a left (often save or run) button and a right cancel button in the specified row, and
        returns them.
        """
        frame = ttk.Frame(self.left)
        rb = RSButton(frame, gui.key.g_cancel, self.cancel)
        lb = LSButton(frame, key, cmd)

        frame.grid(row=row, column=0, columnspan=5, sticky=tk.E)
        rb.place(0)
        lb.place(0)

        return lb, rb


class AnalysisSettings(NWindow):
    """
    Settings window with input fields, check buttons and a save and cancel button.
    """
    def __init__(self, master):
        NWindow.__init__(self, master)
        # title
        self. title = TLabel(self.left, gui.key.as_title)
        # genes
        self.prom = IntInput(self.left, gui.key.as_prom, cfg.user.as_prom, neg=True)
        self.adop = TCheckbutton(self.left, gui.key.ar_adop, cfg.run.a_adop)
        self.pwm = IntInput(self.left, gui.key.as_pwm, cfg.user.as_pwm, pos=True)
        # TFBS
        self.tfbs = IntInput(self.left, gui.key.as_rnd, cfg.user.as_rnd, pos=True)
        # names
        self.int = TextInput(self.left, gui.key.as_int, cfg.user.as_int)
        self.int_short = TextInput(self.left, gui.key.as_int_short, cfg.user.as_int_short)
        self.non_int = TextInput(self.left, gui.key.as_not_int, cfg.user.as_not_int)
        self.non_int_short = TextInput(self.left, gui.key.as_not_int_short, cfg.user.as_not_int_short)
        # run flags
        self.open = TCheckbutton(self.left, gui.key.ar_open, cfg.run.a_open)
        # save and cancel button
        self.lb, self.rb = self.buttons(gui.key.g_save, self.save, 14)
        # place the widgets, active mouse-over and centre the window
        self.place()
        self.mouse_over()
        self.update()
        self.centre()
        # make the main window unresponsive until this window is closed
        master.wait_window(self)

    def place(self):
        """
        Place the widgets in the settings window.
        """
        self.title.place(0)
        self.prom.show(1)
        self.seperator(2)
        self.tfbs.show(3)
        self.pwm.show(4)
        self.seperator(5)
        self.int.show(6)
        self.int_short.show(7)
        self.non_int.show(8)
        self.non_int_short.show(9)
        self.seperator(10)
        self.adop.place(11)
        self.open.place(12)
        self.seperator(13)

    def mouse_over(self):
        """
        Activates mouse-over behaviour for each widget.
        """
        self.title.mouse_over(self.desc)
        self.prom.mouse_over(self.desc)
        self.adop.mouse_over(self.desc)
        self.tfbs.mouse_over(self.desc)
        self.pwm.mouse_over(self.desc)
        self.int.mouse_over(self.desc)
        self.int_short.mouse_over(self.desc)
        self.non_int_short.mouse_over(self.desc)
        self.non_int.mouse_over(self.desc)
        self.open.mouse_over(self.desc)
        self.lb.mouse_over(self.desc)
        self.rb.mouse_over(self.desc)

    def validate(self):
        """
        Checks if the input is valid.
        :return: True if valid, False otherwise
        """
        # the default relative promoter start needs to be a number
        if not self.prom.validate():
            return False

        # the number of random mutation needs to be a positive number
        if not self.tfbs.validate():
            return False

        # the number of required PWM sequences must be a positive number
        if not self.pwm.validate():
            return False

        # genes of interest fields
        if not self.int.validate():
            return False

        if not self.int_short.validate():
            return False

        if not self.non_int.validate():
            return False

        if not self.non_int_short.validate():
            return False

        return True

    def save(self):
        """
        Save the settings and close the window.
        """
        try:
            # don't save and don't close the window if the input isn't valid
            if not self.validate():
                return

            cfg.user.as_prom = int(self.prom.var.get())
            cfg.run.a_adop = self.adop.var.get()
            cfg.user.as_rnd = int(self.tfbs.var.get())
            cfg.run.a_open = self.open.var.get()
            cfg.user.as_pwm = int(self.pwm.var.get())
            cfg.user.as_int_short = self.int_short.var.get()
            cfg.user.as_int = self.int.var.get()
            cfg.user.as_not_int_short = self.non_int_short.var.get()
            cfg.user.as_not_int = self.non_int.var.get()

            cfg.save_config()
            self.destroy()
        except Exception:
            messagebox.showerror('Advanced Mutation Analysis Settings Error',
                                 'An unhandled exception occurred while saving the advanced mutation '
                                 'analysis settings.\n\n'
                                 '{0}'.format(traceback.format_exc()))


class NgsSettings(NWindow):
    """
    Settings window with entry fields, check buttons and a save and cancel button.
    """
    def __init__(self, master):
        NWindow.__init__(self, master)
        # title
        self.title = TLabel(self.left, gui.key.ns_title)
        # executables
        self.bwa = ExeInput(self.left, gui.key.ne_bwa, cfg.user.ne_bwa, 'bwa')
        self.sam = ExeInput(self.left, gui.key.ne_sam, cfg.user.ne_samtools, 'samtools')
        self.var = ExeInput(self.left, gui.key.ne_var, cfg.user.ne_varscan, 'varscan')
        # cut-offs
        self.pval = PercInput(self.left, gui.key.ns_pval, cfg.user.ns_pv)
        self.qual = IntInput(self.left, gui.key.ns_qual, cfg.user.ns_qcut, pos=True)
        # run flags
        self.clean = TCheckbutton(self.left, gui.key.nr_clean, cfg.run.n_clean)
        self.open = TCheckbutton(self.left, gui.key.nr_open, cfg.run.n_open)
        # save and cancel button
        self.lb, self.rb = self.buttons(gui.key.g_save, self.save, 11)
        # place the widgets, active mouse-over and centre the window
        self.place()
        self.mouse_over()
        self.update()
        self.centre()
        # make the main window unresponsive until this window is closed
        master.wait_window(self)

    def place(self):
        """
        Places all widgets in the window.
        """
        self.title.place(0)
        self.bwa.show(1)
        self.sam.show(2)
        self.var.show(3)
        self.seperator(4)
        self.qual.show(5)
        self.pval.show(6)
        self.seperator(7)
        self.clean.place(8)
        self.open.place(9)
        self.seperator(10)

    def mouse_over(self):
        """
        Activates mouse-over behaviour for all widgets.
        """
        self.title.mouse_over(self.desc)
        self.bwa.mouse_over(self.desc)
        self.sam.mouse_over(self.desc)
        self.var.mouse_over(self.desc)
        self.pval.mouse_over(self.desc)
        self.qual.mouse_over(self.desc)
        self.open.mouse_over(self.desc)
        self.clean.mouse_over(self.desc)
        self.lb.mouse_over(self.desc)
        self.rb.mouse_over(self.desc)

    def validate(self):
        """
        Checks if the input is valid.
        :return: True if valid, False otherwise
        """
        # the BWA path needs to point to a valid executable
        if not self.bwa.validate():
            return False

        # the SAMTools path needs to point to a valid executable
        if not self.sam.validate():
            return False

        # the VarScan path needs to point to a valid executable
        if not self.var.validate():
            return False

        # the p-value needs to be a number between 0 and 1
        if not self.pval.validate():
            return False

        # the mapping quality needs to be a positive number
        if not self.qual.validate():
            return False

        return True

    def save(self):
        """
        Save the settings and close the window.
        """
        try:
            # don't save and don't close the window if the input isn't valid
            if not self.validate():
                return

            cfg.user.ne_bwa = adjust_file_path(self.bwa.var.get())
            cfg.user.ne_samtools = adjust_file_path(self.sam.var.get())
            cfg.user.ne_varscan = adjust_file_path(self.var.var.get())

            cfg.user.ns_pv = float(self.pval.var.get())
            cfg.user.ns_qcut = int(self.qual.var.get())

            cfg.run.n_clean = self.clean.var.get()
            cfg.run.n_open = self.open.var.get()

            cfg.save_config()
            self.destroy()
        except Exception:
            messagebox.showerror('Advanced NGS Settings Error',
                                 'An unhandled exception occurred while saving the advanced NGS pipeline settings.\n\n'
                                 '{0}'.format(traceback.format_exc()))


class UniProtDatabaseConverter(NWindow):
    """
    File converter window for converting UniProt database text files.
    """
    def __init__(self, master):
        NWindow.__init__(self, master)
        # title
        self.title = TLabel(self.left, gui.key.pu_title)
        # input and output fields
        self.out_path = SaveFileInput(self.left, gui.key.pu_odir, cfg.user.pu_out)
        self.input = ReadFileInput(self.left, gui.key.pu_in, cfg.user.pu_in)
        # flags
        self.open = TCheckbutton(self.left, gui.key.pu_open, cfg.run.p_open)
        # run and cancel buttons
        self.lb, self.rb = self.buttons(gui.key.pu_run, self.run, 6)
        # place the widgets, active mouse-over and centre the window
        self.place()
        self.mouse_over()
        self.update()
        self.centre()
        # make the main window unresponsive until this window is closed
        master.wait_window(self)

    def place(self):
        """
        Places all widgets in the window.
        """
        self.title.place(0)
        self.out_path.show(1)
        self.input.show(2)
        self.seperator(3)
        self.open.place(4)
        self.seperator(5)

    def mouse_over(self):
        """
        Activates mouse-over behaviour for all widgets.
        """
        self.title.mouse_over(self.desc)
        self.out_path.mouse_over(self.desc)
        self.input.mouse_over(self.desc)
        self.lb.mouse_over(self.desc)
        self.rb.mouse_over(self.desc)
        self.open.mouse_over(self.desc)

    def validate(self):
        """
        Checks if the input is valid.
        :return: True if valid, False otherwise
        """
        # the output file needs to exist
        if not self.out_path.validate(exists=False):
            return False

        # the input directory needs to exist
        if not self.input.validate():
            return False

        return True

    def save(self):
        """
        Save the input and output paths to the user configuration.
        """
        cfg.user.pu_out = self.out_path.var.get()
        cfg.user.pu_in = self.input.var.get()

        cfg.save_uniprot()

    def run(self):
        """
        Convert a UniProt protein text file to a protein domain .tsv file.
        """
        try:
            # do nothing if the input is not valid
            if not self.validate():
                return

            self.save()

            in_path = adjust_dir_path(self.input.var.get())
            out_path = adjust_file_path(self.out_path.var.get())

            UniProtConverter(in_path, out_path)

            # open the directory of the output file
            if self.open.var.get():
                self.out_path.open()
        except Exception:
            messagebox.showerror('UniProt Converter Error',
                                 'An unhandled exception occurred while converting the UniProt file.\n\n'
                                 '{0}'.format(traceback.format_exc()))


class PatricARConverter(NWindow):
    """
    File converter window for converting UniProt database text files.
    """
    def __init__(self, master):
        NWindow.__init__(self, master)
        # title
        self.title = TLabel(self.left, gui.key.pp_title)
        # input and output fields
        self.out_path = SaveFileInput(self.left, gui.key.pp_odir, cfg.user.pp_out)
        self.input = ReadFileInput(self.left, gui.key.pp_in, cfg.user.pp_in)
        # flags
        self.open = TCheckbutton(self.left, gui.key.pp_open, cfg.run.p_open)
        # run and cancel buttons
        self.lb, self.rb = self.buttons(gui.key.pp_run, self.run, 6)
        # place the widgets, active mouse-over and centre the window
        self.place()
        self.mouse_over()
        self.update()
        self.centre()
        # make the main window unresponsive until this window is closed
        master.wait_window(self)

    def place(self):
        """
        Places all widgets in the window.
        """
        self.title.place(0)
        self.out_path.show(1)
        self.input.show(2)
        self.seperator(3)
        self.open.place(4)
        self.seperator(5)

    def mouse_over(self):
        """
        Activates mouse-over behaviour for all widgets.
        """
        self.title.mouse_over(self.desc)
        self.out_path.mouse_over(self.desc)
        self.input.mouse_over(self.desc)
        self.lb.mouse_over(self.desc)
        self.rb.mouse_over(self.desc)
        self.open.mouse_over(self.desc)

    def validate(self):
        """
        Checks if the input is valid.
        :return: True if valid, False otherwise
        """
        # the output file needs to exist
        if not self.out_path.validate(exists=False):
            return False

        # the input directory needs to exist
        if not self.input.validate():
            return False

        return True

    def save(self):
        """
        Save the input and output paths to the user configuration.
        """
        cfg.user.pu_out = self.out_path.var.get()
        cfg.user.pu_in = self.input.var.get()

        cfg.save_patric()

    def run(self):
        """
        Convert PATRIC antibiotic resistance file to a .tsv file.
        """
        try:
            # do nothing if the input is not valid
            if not self.validate():
                return

            in_path = adjust_dir_path(self.input.var.get())
            out_path = adjust_file_path(self.out_path.var.get())

            PatricConverter(in_path, out_path)

            # open the directory of the output file
            if self.open.var.get():
                self.out_path.open()
        except Exception:
            messagebox.showerror('PATRIC Converter Error',
                                 'An unhandled exception occurred while converting the PATRIC '
                                 'antibiotic resistance file.\n\n'
                                 '{0}'.format(traceback.format_exc()))


class VcfMerger(NWindow):
    """
    File converter window for merging .vcf files with mutations.
    """
    def __init__(self, master):
        NWindow.__init__(self, master)
        # title
        self.title = TLabel(self.left, gui.key.pm_title)
        # input and output fields
        self.out_path = SaveFileInput(self.left, gui.key.pm_odir, cfg.user.pm_out)
        self.input = DirInput(self.left, gui.key.pm_in, cfg.user.pm_in)
        # flags
        self.open = TCheckbutton(self.left, gui.key.pm_open, cfg.run.p_open)
        # run and cancel buttons
        self.lb, self.rb = self.buttons(gui.key.pm_run, self.run, 6)
        # place the widgets, active mouse-over and centre the window
        self.place()
        self.mouse_over()
        self.update()
        self.centre()
        # make the main window unresponsive until this window is closed
        master.wait_window(self)

    def place(self):
        """
        Places all widgets in the window.
        """
        self.title.place(0)
        self.out_path.show(1)
        self.input.show(2)
        self.seperator(3)
        self.open.place(4)
        self.seperator(5)

    def mouse_over(self):
        """
        Activates mouse-over behaviour for all widgets.
        """
        self.title.mouse_over(self.desc)
        self.out_path.mouse_over(self.desc)
        self.input.mouse_over(self.desc)
        self.lb.mouse_over(self.desc)
        self.rb.mouse_over(self.desc)
        self.open.mouse_over(self.desc)

    def validate(self):
        """
        Checks if the input is valid.
        :return: True if valid, False otherwise
        """
        # the output file needs to exist
        if not self.out_path.validate(exists=False):
            return False

        # the input directory needs to exist
        if not self.input.validate():
            return False

        return True

    def save(self):
        """
        Save the input and output paths to the user configuration.
        """
        cfg.user.pm_out = self.out_path.var.get()
        cfg.user.pm_in = self.input.var.get()

        cfg.save_merger()

    def run(self):
        """
        Merge mutation .vcf files into a single .tsv file.
        """
        try:
            # do nothing if the input is not valid
            if not self.validate():
                return

            self.save()

            in_path = adjust_dir_path(self.input.var.get())
            out_path = adjust_file_path(self.out_path.var.get())

            m = MutationVcfMerger(in_path, out_path)

            # open the directory of the output file if specified and no error occurred
            if self.open.var.get() and not m.empty and not m.no_files and not m.write_error:
                self.out_path.open()
        except Exception:
            messagebox.showerror('Mutation VCF Merger Error',
                                 'An unhandled exception occurred while merging mutation VCF files.\n\n'
                                 '{0}'.format(traceback.format_exc()))


class RegulonConverter(NWindow):
    """
    File converter window for merging .vcf files with mutations.
    """
    def __init__(self, master):
        NWindow.__init__(self, master)
        # title
        self.title = TLabel(self.left, gui.key.pr_title)
        # input and output fields
        self.out_path = DirInput(self.left, gui.key.pr_odir, cfg.user.pr_out)

        self.gene = ReadFileInput(self.left, gui.key.pr_gene, cfg.user.pr_gene)
        self.operon = ReadFileInput(self.left, gui.key.pr_op, cfg.user.pr_operon)
        self.tu = ReadFileInput(self.left, gui.key.pr_tu, cfg.user.pr_tu)

        self.tu_reg = ReadFileInput(self.left, gui.key.pr_tu_reg, cfg.user.pr_tftu)
        self.gene_reg = ReadFileInput(self.left, gui.key.pr_gene_reg, cfg.user.pr_tg)
        self.tf_reg = ReadFileInput(self.left, gui.key.pr_tf_reg, cfg.user.pr_tftf)
        self.operon_reg = ReadFileInput(self.left, gui.key.pr_operon, cfg.user.pr_to)

        self.msa = ReadFileInput(self.left, gui.key.pr_msa, cfg.user.pr_msa)
        self.tfbs = ReadFileInput(self.left, gui.key.pr_tfbs, cfg.user.pr_tfbs)
        # status window
        self.status = Status(self.right, self, self.desc, 2)
        # flags
        self.skip = TCheckbutton(self.left, gui.key.pr_skip, False)
        self.open = TCheckbutton(self.left, gui.key.pr_open, cfg.run.p_open)
        # run and cancel buttons
        self.lb, self.rb = self.buttons(gui.key.pr_run, self.run, 26)
        # place the widgets, active mouse-over and centre the window
        self.place()
        self.mouse_over()
        self.update()
        self.centre()
        # make the main window unresponsive until this window is closed
        master.wait_window(self)

    def place(self):
        """
        Places all widgets in the window.
        """
        self.title.place(0)
        self.out_path.show(1)
        self.seperator(2)
        self.gene.show(3)
        self.operon.show(5)
        self.tu.show(6)
        self.seperator(7)
        self.gene_reg.show(8)
        self.tf_reg.show(9)
        self.operon_reg.show(10)
        self.tu_reg.show(11)
        self.seperator(12)
        self.tfbs.show(18)
        self.msa.show(20)
        self.seperator(22)
        self.skip.place(23)
        self.open.place(24)
        self.seperator(25)

    def mouse_over(self):
        """
        Activates mouse-over behaviour for all widgets.
        """
        self.title.mouse_over(self.desc)
        self.out_path.mouse_over(self.desc)
        self.gene.mouse_over(self.desc)
        self.operon.mouse_over(self.desc)
        self.tu.mouse_over(self.desc)
        self.tu_reg.mouse_over(self.desc)
        self.operon_reg.mouse_over(self.desc)
        self.gene_reg.mouse_over(self.desc)
        self.tf_reg.mouse_over(self.desc)
        self.tfbs.mouse_over(self.desc)
        self.msa.mouse_over(self.desc)
        self.skip.mouse_over(self.desc)
        self.lb.mouse_over(self.desc)
        self.rb.mouse_over(self.desc)
        self.open.mouse_over(self.desc)

    def save(self):
        """
        Save the current input paths to the user configuration.
        """
        cfg.user.pr_out = self.out_path.var.get()
        cfg.user.pr_gene = self.gene.var.get()
        cfg.user.pr_operon = self.operon.var.get()
        cfg.user.pr_tu = self.tu.var.get()
        cfg.user.pr_tftu = self.tu_reg.var.get()
        cfg.user.pr_tftf = self.tf_reg.var.get()
        cfg.user.pr_tg = self.gene_reg.var.get()
        cfg.user.pr_to = self.operon_reg.var.get()
        cfg.user.pr_tfbs = self.tfbs.var.get()
        cfg.user.pr_msa = self.msa.var.get()

        cfg.save_regulon()

    def validate(self):
        """
        Checks if the input is valid.
        :return: True if valid, False otherwise
        """
        # the output directory needs to exist
        if not self.setup():
            return False
        if not self.out_path.validate():
            return False

        # the input file needs to exist
        if not self.gene.validate():
            return False

        # the input file needs to exist
        if not self.gene_reg.validate():
            return False

        # the input file needs to exist
        if not self.operon_reg.validate():
            return False

        # the input file needs to exist
        if not self.tf_reg.validate():
            return False

        # the input file needs to exist
        if not self.tfbs.validate():
            return False

        # the input file needs to exist
        if not self.msa.validate():
            return False

        # the input file needs to exist
        if not self.operon.validate():
            return False

        # the input file needs to exist
        if not self.tu.validate():
            return False

        # the input file needs to exist
        if not self.tu_reg.validate():
            return False

        return True

    def run(self):
        """
        Merge mutation .vcf files into a single .tsv file.
        """
        try:
            # do nothing if the input is not valid
            if not self.validate():
                return

            self.save()

            # set up the converter arguments
            RegulonDBConverter(out_dir=self.out_path.var.get(),
                               gene=self.gene.var.get(),
                               operons=self.operon.var.get(),
                               tus=self.tu.var.get(),
                               operon_reg=self.operon_reg.var.get(),
                               gene_reg=self.gene_reg.var.get(),
                               tf_reg=self.tf_reg.var.get(),
                               tu_reg=self.tu_reg.var.get(),
                               tfbs=self.tfbs.var.get(),
                               msa=self.msa.var.get(),
                               skip=self.skip.var.get(),
                               status=self.status)

            # open the directory of the output file
            if self.open.var.get():
                self.out_path.open()
        except Exception:
            messagebox.showerror('RegulonDB Converter Error',
                                 'An unhandled exception occurred while converting RegulonDB files.\n\n'
                                 '{0}'.format(traceback.format_exc()))


class AnalysisSection:
    """
    Creates and sets up a section for mutation analysis with input fields and run button.
    """
    def __init__(self, master, desc, status, row):
        self.master = master
        self.title = TLabel(master, gui.key.ah_title)

        self.out_dir = DirInput(master, gui.key.ao_dir, cfg.user.ao_dir)

        self.genes = ReadFileInput(master, gui.key.ai_gene, cfg.user.ai_gene)
        self.mut = ReadFileInput(master, gui.key.ai_mut, cfg.user.ai_mut)

        self.int = ReadFileInput(master, gui.key.ai_int, cfg.user.ai_int)
        self.int_sel = Selector(master, gui.key.ah_int, cfg.run.a_int, [self.int])

        self.sm = ReadFileInput(master, gui.key.ai_sm, cfg.user.ai_sm)
        self.sm_sel = Selector(master, gui.key.ah_coding, cfg.run.a_coding, [self.sm])

        self.prot = ReadFileInput(master, gui.key.ai_prot, cfg.user.ai_prot)
        self.prot_sel = Selector(master, gui.key.ah_prot, cfg.run.a_prot, [self.prot])

        self.reg = ReadFileInput(master, gui.key.ai_reg, cfg.user.ai_reg)
        self.reg_sel = Selector(master, gui.key.ah_reg, cfg.run.a_reg, [self.reg])

        self.tfbs = ReadFileInput(master, gui.key.ai_tfbs, cfg.user.ai_tfbs)
        self.msa = ReadFileInput(master, gui.key.ai_msa, cfg.user.ai_msa)
        self.tfbs_sel = Selector(master, gui.key.ah_tfbs, cfg.run.a_tfbs, [self.tfbs, self.msa])

        self.run_button = RButton(master, gui.key.ar_run, lambda: self.run(status))

        self.place(row)
        self.mouse_over(desc)

    def place(self, row):
        """
        Places the analysis section widgets in the main application window.
        :param row: start row of the section
        """
        self.title.place(row)
        # output directory
        self.out_dir.show(row + 1)
        # separator between required output and input
        sep = ttk.Separator(self.master, orient=tk.HORIZONTAL)
        sep.grid(row=row + 2, column=1, columnspan=4, pady=gui.sp_pady, padx=gui.sp_padx, sticky=tk.EW)
        # required genes and mutation input
        self.genes.show(row + 3)
        self.mut.show(row + 4)
        # optional sections_verbose with their separators, check buttons and input fields
        self.int_sel.place(row + 5)
        self.sm_sel.place(row + 7)
        self.prot_sel.place(row + 9)
        self.reg_sel.place(row + 11)
        self.tfbs_sel.place(row + 13)
        # separator between the input fields and the run button
        sep = ttk.Separator(self.master, orient=tk.HORIZONTAL)
        sep.grid(row=row + 16, column=1, columnspan=4, pady=gui.sp_pady, padx=gui.sp_padx, sticky=tk.EW)
        # run button
        self.run_button.place(row + 17)

    def mouse_over(self, desc):
        """
        Activates mouse-over behaviour to each widget of the analysis section.
        :param desc: description box
        """
        self.title.mouse_over(desc)
        self.out_dir.mouse_over(desc)
        self.genes.mouse_over(desc)
        self.mut.mouse_over(desc)

        self.int_sel.mouse_over(desc)
        self.sm_sel.mouse_over(desc)
        self.prot_sel.mouse_over(desc)
        self.reg_sel.mouse_over(desc)
        self.tfbs_sel.mouse_over(desc)

        self.run_button.mouse_over(desc)

    def save(self):
        """
        Updates the user configuration with values from the entry fields.
        """
        cfg.user.ao_dir = adjust_dir_path(self.out_dir.var.get())
        cfg.user.ai_gene = adjust_file_path(self.genes.var.get())
        cfg.user.ai_mut = adjust_file_path(self.mut.var.get())
        cfg.user.ai_int = adjust_file_path(self.int.var.get())
        cfg.user.ai_reg = adjust_file_path(self.reg.var.get())
        cfg.user.ai_sm = adjust_file_path(self.sm.var.get())
        cfg.user.ai_prot = adjust_file_path(self.prot.var.get())
        cfg.user.ai_tfbs = adjust_file_path(self.tfbs.var.get())
        cfg.user.ai_msa = adjust_file_path(self.msa.var.get())

        cfg.run.a_int = self.int_sel.var.get()
        cfg.run.a_reg = self.reg_sel.var.get()
        cfg.run.a_coding = self.sm_sel.var.get()
        cfg.run.a_prot = self.prot_sel.var.get()
        cfg.run.a_tfbs = self.tfbs_sel.var.get()

    def update(self):
        """
        Updates the entry fields with the current values from the user configuration.
        """
        self.out_dir.var.set(cfg.user.ao_dir)
        self.genes.var.set(cfg.user.ai_gene)
        self.mut.var.set(cfg.user.ai_mut)
        self.int.var.set(cfg.user.ai_int)
        self.reg.var.set(cfg.user.ai_reg)
        self.sm.var.set(cfg.user.ai_sm)
        self.prot.var.set(cfg.user.ai_prot)
        self.tfbs.var.set(cfg.user.ai_tfbs)
        self.msa.var.set(cfg.user.ai_msa)
        
        self.int_sel.var.set(cfg.run.a_int)
        self.reg_sel.var.set(cfg.run.a_reg)
        self.sm_sel.var.set(cfg.run.a_coding)
        self.prot_sel.var.set(cfg.run.a_prot)
        self.tfbs_sel.var.set(cfg.run.a_tfbs)

        self.master.update()

    def validate(self):
        """
        Checks if the analysis input is valid and updates the configuration.
        :return: True if the input is valid and the analysis can be run, False otherwise
        """
        # make sure the output directory exists or is created. if that is not possible, give an error message
        # and don'r run the analysis
        if not self.out_dir.setup():
            return False
        if not self.out_dir.clean():
            return False

        # the gene.tsv file is required for the analysis. if it is not given or not accessible for reading,
        # give an error message and don't run the analysis
        if not self.genes.validate():
            return False

        # the mutation.tsv file is required for the analysis. if it is not given or not accessible for reading,
        # give an error message and don't run the analysis
        if not self.mut.validate():
            return False

        # the substitution matrix file is required for the coding region analysis. if it is not given or not
        # accessible for reading, don't run the analysis
        if self.sm_sel.var.get() and not self.sm.validate():
            return False

        # the protein domain .tsv file is required for adding protein domain information. if it is not given or not
        # accessible for reading, don't run the protein domain parsing
        if self.prot_sel.var.get() and not self.prot.validate():
            return False

        # the interest.tsv file is required for adding genes of interests. if it is not given or not
        # accessible for reading, don't run the analysis
        if self.int_sel.var.get() and not self.int.validate():
            return False

        # the regulation.tsv file is required for adding regulatory information. if it is not given or not accessible
        # for reading, add a warning to the log and don't run the analysis
        if self.reg_sel.var.get() and not self.reg.validate():
            return False

        if self.tfbs_sel.var.get():
            # the TFBS.tsv file is required for analysing TF binding sites. if it is not given or not accessible
            # for reading, don't run the analysis
            if not self.tfbs.validate():
                return False
            # the TF MSA .tsv file is required for analysing TF binding sites. if it is not given or not accessible
            # for reading, don't run the analysis
            if not self.msa.validate():
                return False

        # update the configuration
        self.save()
        # update the output fields
        cfg.update_output(cfg.user.ao_dir, cfg.user.no_dir)

        return True

    def run(self, status):
        """
        Executes all analysis steps that can be performed based on the input.
        :type status: Status
        """
        try:
            # clear the previous status messages
            status.clear()

            if not out_warn(self.out_dir.var.get()):
                return

            # header status line
            status.name('Mutation Analysis')

            # start the runtime clock
            start = time()
            stop = False

            status.running(0, 'pre-processing input')

            # reset the log
            cfg.a_log = ALog()

            # make sure the input is valid
            if not self.validate():
                status.fail(0, 'pre-processing input')
                stop = True

            # reset sub-categories
            cfg.cat = set()

            # build database
            db = Database()

            if not stop:
                worked = db.build_database()

                if len(cfg.cat) > 10:
                    messagebox.showerror('Mutation Analysis Error',
                                         'The number of sub-categories in the "genes of interest" file must be '
                                         'smaller or equal 10.')

                if worked:
                    status.time(0, 'pre-processing input', time() - start)
                else:
                    status.fail(0, 'pre-processing input')
                    stop = True

            # set up the statistics
            stats = Stats()

            # analyse coding region
            if cfg.run.a_coding and not stop:
                start2 = time()
                status.running(1, 'analysing coding regions')
                coding_scores = db.analyse_coding_region()

                if db.coding:
                    stats.set_coding_scores(coding_scores)
                    status.time(1, 'analysing coding regions', time() - start2)
                else:
                    status.fail(1, 'analysing coding regions')
            else:
                status.skip(1, 'analysing coding regions')

            # analyse transcription factor binding sites
            if cfg.run.a_tfbs and not stop:
                start2 = time()
                status.running(2, 'analysing TFBS')
                tfbs_scores = db.analyse_tfbs()
                stats.set_tfbs_scores(tfbs_scores)
                status.time(2, 'analysing TFBS', time() - start2)
            else:
                status.skip(2, 'analysing TFBS')

            # compute general statistics
            if not stop:
                start2 = time()
                status.running(3, 'computing general statistics')
                stats.compute(db)
                status.time(3, 'computing general statistics', time() - start2)
            else:
                status.skip(3, 'computing general statistics')

            # generate plot
            if cfg.run.a_int and not stop:
                start2 = time()
                status.running(4, 'generating plots')
                stats.plot()
                status.time(4, 'generating plots', time() - start2)
            else:
                status.skip(4, 'generating plots')

            # write statistics
            if not stop:
                start2 = time()
                status.running(5, 'writing statistics')
                worked = stats.write()

                if worked:
                    status.time(5, 'writing statistics', time() - start2)
                else:
                    status.fail(5, 'writing statistics')
            else:
                status.skip(5, 'writing statistics')

            # build GRN
            if cfg.run.a_reg and not stop:
                start2 = time()
                status.running(6, 'generating GRN')
                db.build_grn()
                status.time(6, 'generating GRN', time() - start2)
            else:
                status.skip(6, 'generating GRN')

            # write known mutation list
            if cfg.run.a_prot and cfg.run.a_int and not stop:
                start2 = time()
                status.running(7, 'writing known mutations')
                worked = db.write_mutation_list()

                if worked:
                    status.time(7, 'writing known mutations', time() - start2)
                else:
                    status.fail(7, 'writing known mutations')
            else:
                status.skip(7, 'writing known mutations')

            # write gene mutation list
            if cfg.run.a_int and not stop:
                start2 = time()
                status.running(8, 'writing genes of interest table')
                worked = db.write_gene_list()

                if worked:
                    status.time(8, 'writing genes of interest table', time() - start2)
                else:
                    status.fail(8, 'writing genes of interest table')
            else:
                status.skip(8, 'writing genes of interest table')

            # write database
            if not stop:
                start2 = time()
                status.running(9, 'writing database')
                worked = db.write_database()

                if worked:
                    status.time(9, 'writing database', time() - start2)
                else:
                    status.fail(9, 'writing database')
            else:
                status.skip(9, 'writing database')

            # write log
            start2 = time()
            status.running(10, 'writing log')
            worked = cfg.a_log.write(cfg.op.a_log)

            if worked:
                status.time(10, 'writing log', time() - start2)
            else:
                status.fail(10, 'writing log')

            # write the total running time
            status.total_time(time() - start)

            # open the output directory if selected
            if cfg.run.a_open:
                self.out_dir.open()
        except Exception:
            messagebox.showerror('Mutation Analysis Error',
                                 'An unhandled exception occurred during the mutation analysis.\n\n'
                                 '{0}'.format(traceback.format_exc()))


class NgsSection:
    """
    Creates and sets up a section for the NGS pipeline with input fields and run button.
    """
    def __init__(self, master, desc, status, row):
        self.master = master
        self.title = TLabel(master, gui.key.nh_title)

        self.out_dir = DirInput(master, gui.key.no_dir, cfg.user.no_dir)

        self.ref = ReadFileInput(master, gui.key.ni_ref, cfg.user.ni_ref)
        self.reads = DirInput(master, gui.key.ni_reads, cfg.user.ni_reads)

        self.run_button = RButton(master, gui.key.nr_run, lambda: self.run(status))

        self.mouse_over(desc)
        self.place(row)

    def place(self, row):
        """
        Places the NGS section widgets in the main application window.
        :param row: start row of the section
        """
        self.title.place(row)
        self.out_dir.show(row + 1)

        sep = ttk.Separator(self.master, orient=tk.HORIZONTAL)
        sep.grid(row=row + 2, column=1, columnspan=4, pady=gui.sp_pady, padx=gui.sp_padx, sticky=tk.EW)

        self.ref.show(row + 3)
        self.reads.show(row + 4)

        sep = ttk.Separator(self.master, orient=tk.HORIZONTAL)
        sep.grid(row=row + 5, column=1, columnspan=4, pady=gui.sp_pady, padx=gui.sp_padx, sticky=tk.EW)

        self.run_button.place(row + 6)

    def mouse_over(self, desc):
        """
        Activates mouse-over behaviour to each widget of the NGS section.
        :param desc: description box
        """
        self.title.mouse_over(desc)
        self.out_dir.mouse_over(desc)
        self.ref.mouse_over(desc)
        self.reads.mouse_over(desc)
        self.run_button.mouse_over(desc)

    def save(self):
        """
        Updates the user configuration with values from the entry fields.
        """
        cfg.user.ni_reads = adjust_dir_path(self.reads.var.get())
        cfg.user.ni_ref = adjust_file_path(self.ref.var.get())
        cfg.user.no_dir = adjust_dir_path(self.out_dir.var.get())

    def update(self):
        """
        Updates the entry fields with the current values from the user configuration.
        """
        self.reads.var.set(cfg.user.ni_reads)
        self.ref.var.set(cfg.user.ni_ref)
        self.out_dir.var.set(cfg.user.no_dir)
        self.master.update()

    def validate(self):
        """
        Checks if the analysis input is valid and updates the configuration.
        :return: True if the input is valid and the analysis can be run, False otherwise
        """
        # make sure the output directory exists or is created. if that is not possible, give an error message
        # and don'r run the NGS pipeline
        if not self.out_dir.setup():
            return False
        if not self.out_dir.clean():
            return False

        # the reference genome .fasta file is required for the NGS pipeline. if it is not given or not accessible for
        # reading, give an error message and don't run the NGS pipeline
        if not self.ref.validate():
            return False

        # the reads directory is required for the NGS pipeline. if it is not a valid directory, give an error message
        # and don't rn the NGS pipeline
        if not self.reads.validate():
            return False

        # update the configuration
        self.save()
        # update the output fields
        cfg.update_output(cfg.user.ao_dir, cfg.user.no_dir)
    
        return True

    def run(self, status):
        """
        Executes the NGS pipeline.
        :type status: Status
        """
        try:
            # clear the previous status messages
            status.clear()

            # warn the user that old output is going to be deleted
            if not out_warn(self.out_dir.var.get()):
                return

            # start the run time clock
            start = time()

            status.name('NGS Pipeline')
            status.running(0, 'pre-processing input')

            # stop if parts of the input are invalid
            if not self.validate():
                status.fail(0, 'pre-processing input')
                return

            # initialise the pipeline based on the input and the application configuration
            pipeline = NGSPipeline()

            # stop if the reference genome is not given as a .fasta file
            if pipeline.ref_not_fasta:
                messagebox.showerror('Error', 'The reference genome file is not a .fasta file. '
                                              '({0})'.format(cfg.user.ni_ref))
                status.fail(0, 'pre-processing input')
                return

            # stop if there are no .fastq files in the reads directory
            if pipeline.reads_no_fastq:
                messagebox.showerror('Error', 'There are no valid .fastq files in the reads directory. '
                                              '({0})'.format(cfg.user.ni_reads))
                status.fail(0, 'pre-processing input')
                return

            # stop if there is not a single strain with the required 2 .fastq files
            if pipeline.reads_no_strains:
                messagebox.showerror('Error', 'There are no strains with two .fastq files in the reads directory. '
                                              '({0})'.format(cfg.user.ni_reads))
                status.fail(0, 'pre-processing input')
                return

            # stop if the read file types are inconsistent
            if pipeline.inconsistent_suffixes:
                msg = 'The file types in the reads directory are inconsistent:\n- {0}'
                messagebox.showerror('Error', msg.format('\n- '.join(sorted(pipeline.inconsistent_suffixes))))
                status.fail(0, 'pre-processing input')
                return

            # pre-processing is finished
            status.time(0, 'pre-processing input', time() - start)

            # execute the actual pipeline
            pipeline.run_pipeline(status)

            # update the GUI with the total running time
            status.total_time(time() - start)

            # open the output directory if selected
            if cfg.run.n_open:
                self.out_dir.open()
        except Exception:
            messagebox.showerror('NGS Pipeline Error',
                                 'An unhandled exception occurred during the NGS pipeline.\n\n'
                                 '{0}'.format(traceback.format_exc()))


class GUI:
    def __init__(self, app):
        self.app = app
        # set up the widget styles
        self.style = ttk.Style()
        self.styles()
        # set up the frames
        self.left, self.right = self.setup()
        # set up the sections_verbose
        self.desc = Description(self.right, 0)
        self.desc.configure(height=gui.df_height, width=gui.df_width)
        self.status = Status(self.right, self.app, self.desc, 2)
        self.ana = AnalysisSection(self.left, self.desc, self.status, 9)
        self.ngs = NgsSection(self.left, self.desc, self.status, 0)
        # set up the menu bar
        self.menubar()
        # centre the application window
        self.centre()

    def menubar(self):
        """
        Sets up the menu bar with sub-menus for settings, the file parser and help.
        """
        """ MENU BAR """
        menubar = tk.Menu(self.app)
        # set the menu as the menu bar of the application window
        self.app.configure(menu=menubar)

        """ SETTINGS MENU"""
        # settings menu
        settings = tk.Menu(menubar, tearoff=0)
        # save settings, restore default settings and restore + save default settings
        settings.add_command(label=gui.desc[gui.key.ms_save].label, command=self.save)
        settings.add_command(label=gui.desc[gui.key.ms_def].label, command=self.default)
        settings.add_command(label=gui.desc[gui.key.ms_save_def].label, command=self.save_default)
        settings.add_separator()
        # advanced settings for analysis and NGS
        settings.add_command(label=gui.desc[gui.key.ms_ngs].label, command=lambda: NgsSettings(self.app))
        settings.add_command(label=gui.desc[gui.key.ms_ana].label, command=lambda: AnalysisSettings(self.app))
        # add the settings menu to the menu bar
        menubar.add_cascade(label=gui.desc[gui.key.ms_title].label, menu=settings)

        """ PARSER MENU """
        # parser menu
        parser = tk.Menu(menubar, tearoff=0)
        parser.add_command(label=gui.desc[gui.key.mc_vmerge].label, command=lambda: VcfMerger(self.app))
        parser.add_command(label=gui.desc[gui.key.mc_uconvert].label,
                           command=lambda: UniProtDatabaseConverter(self.app))
        parser.add_command(label=gui.desc[gui.key.mc_patric].label,
                           command=lambda: PatricARConverter(self.app))
        parser.add_command(label=gui.desc[gui.key.mc_regulon].label,
                           command=lambda: RegulonConverter(self.app))
        # add the parser menu to the menu bar
        menubar.add_cascade(label=gui.desc[gui.key.mc_title].label, menu=parser)

        """ HELP MENU """
        # help menu
        help = tk.Menu(menubar, tearoff=0)
        # open manual
        help.add_command(label=gui.desc[gui.key.mh_man].label, command=lambda: self.open(cfg.manual))
        # open installation guide
        help.add_command(label=gui.desc[gui.key.mh_inst].label, command=lambda: self.open(cfg.installation))
        # add the help menu to the menu bar
        menubar.add_cascade(label=gui.desc[gui.key.mh_title].label, menu=help)

        return menubar

    def setup(self):
        """
        Sets up window borders and section frames.
        :return: left frame, right frame
        """
        # set the title
        self.app.title(gui.title)
        # set up the window border frames and the separating frame
        ttk.Frame(self.app).grid(row=0, column=0, sticky=tk.NSEW, ipadx=gui.wl_padx)
        ttk.Frame(self.app).grid(row=0, column=4, sticky=tk.NSEW, ipadx=gui.wr_padx)
        ttk.Frame(self.app).grid(row=0, column=2, sticky=tk.NSEW, ipadx=gui.col_padx)
        # don't show bottom padding on Mac OS
        if cfg.os != cfg.mac:
            ttk.Frame(self.app).grid(row=1, column=0, ipady=gui.wb_pady)
        self.app.rowconfigure(0, weight=1)
        # set up the left and right section frame
        left = ttk.Frame(self.app)
        right = ttk.Frame(self.app)
        left.grid(row=0, column=1, sticky=tk.NSEW)
        right.grid(row=0, column=3, sticky=tk.NSEW)
        # change the expanding behaviour of the right frame
        right.rowconfigure(1, weight=1)
        left.rowconfigure(7, weight=1)

        return left, right

    def centre(self):
        """
        Places the application window in the middle of the screen.
        """
        self.app.update()
        # compute the current height and width
        h = self.app.winfo_reqheight()
        w = self.app.winfo_reqwidth()
        # compute the x and y coordinates of the window using the screen width and height
        x = (self.app.winfo_screenwidth() - w) // 2
        y = (self.app.winfo_screenheight() - h) // 2

        # set the window position (x coordinate, y coordinate) on the screen
        self.app.geometry('{0}x{1}+{2}+{3}'.format(w, h, x, y))
        self.app.update()

    def fonts(self, style, fnt, fg=None):
        """
        Build the font by trying to convert it into a Tkinter font.
        :param style: base_style name
        :param fnt: tuple of font name, and optionally size and weight
        :param fg: text colour
        """
        # no font specified
        if not fnt:
            return

        if not fnt[0]:
            return
        # see if the font is a TKinter font
        try:
            temp = (font.nametofont(fnt[0]), *fnt[1:])
        except Exception:
            temp = fnt

        if fg:
            self.style.configure(style, foreground=fg, font=temp)
        else:
            self.style.configure(style, font=temp)

    def styles(self):
        """
        Sets up the widget styles.
        """
        if cfg.gui.theme in self.style.theme_names():
            self.style.theme_use(cfg.gui.theme)

        self.fonts(gui.dt_style, gui.dt_font, gui.dt_colour)
        self.fonts(gui.st_style, gui.st_font, gui.st_colour)
        self.fonts(gui.sst_style, gui.sst_font, gui.sst_colour)
        self.fonts(gui.fb_style, gui.fb_font)
        self.fonts(gui.el_style, gui.el_font)
        self.fonts(gui.sl_style, gui.sl_font)
        self.fonts(gui.cb_style, gui.cb_font)
        self.fonts(gui.rb_style, gui.rb_font)

    def default(self):
        """
        Restore the default settings.
        """
        cfg.restore_default()
        self.ana.update()
        self.ngs.update()

    def save(self):
        """
        Save the current settings.
        """
        self.ana.save()
        self.ngs.save()
        cfg.save_config()

    def save_default(self):
        """
        First restore default settings and then save them.
        """
        self.default()
        self.save()

    @staticmethod
    def open(file):
        """
        Opens the user manual PDF file.
        """
        try:
            # open command depending on the operating system
            if cfg.os == cfg.linux:
                Popen(["xdg-open", file])
            elif cfg.os == cfg.mac:
                Popen(["open", file])
            elif cfg.os == cfg.windows:
                os.startfile(file)
            elif cfg.os == cfg.cygwin:
                os.startfile(file)
        except Exception:
            messagebox.showerror('Open Error',
                                 'An unhandled exception occurred when opening the following file:\n{0}'
                                 '\n\n{1}'.format(file, traceback.format_exc()))
