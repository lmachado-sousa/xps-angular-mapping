import os
import re
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import matplotlib as mpl
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

import angular_map
import colormaps
import gaussian
import savitzky
from shirley import shirley

class MainWindow(tk.Frame):
    def __init__(self,parent=None):
#Creating frames
        tk.Frame.__init__(self)
        self.pack(expand=True,fill='both')
        self.menu_frame = tk.Frame(self)
        self.menu_frame.pack(side=tk.TOP)
        self.image_frame = tk.Frame(self)
        self.image_frame.pack(side=tk.TOP, expand=True, fill='both')
        self.inter_frame = tk.Frame(self)
        self.inter_frame.pack(expand=True)
        self.info_frame = tk.Frame(self)
        self.info_frame.pack(side=tk.TOP, expand=True, fill='both')
        self.slider_frame = tk.Frame(self)
        self.slider_frame.pack(side=tk.TOP)
        self.progress_frame = tk.Frame(self)
        self.progress_frame.pack(side=tk.TOP)
#Change Tkinter themes using ttk
        self.style = ttk.Style()
        self.style.configure('TButton',relief='flat')
        self.available_themes = self.style.theme_names()
#Creating widgets
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig,master=self.image_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side='bottom', fill=tk.BOTH, expand=True)
        self.label = tk.Label(self.info_frame)
        self.label.pack(side=tk.LEFT)
        self.evlabel = tk.Label(self.info_frame)
        self.evlabel.pack(side=tk.RIGHT)
        self.slider = tk.Scale(self.slider_frame,orient=tk.HORIZONTAL,length=300,
                                    command=self.slide_image,to=0,sliderrelief='flat')
        self.slider.pack(side=tk.BOTTOM)
        self.progress = ttk.Progressbar(self.progress_frame,orient="horizontal",length=300, mode="determinate")
        self.progress.pack()

#Creates main application buttons
        self.openbutton = ttk.Button(self.menu_frame, text = "Open", command=self.openfile).pack(side=tk.LEFT)
        self.prevbutton = ttk.Button(self.inter_frame, text = "<", command=self.previmg, state=tk.DISABLED)
        self.prevbutton.pack(side=tk.LEFT)
        self.nextbutton = ttk.Button(self.inter_frame, text = ">", command=self.nextimg, state=tk.DISABLED)
        self.nextbutton.pack(side=tk.LEFT)
        self.testbutton = ttk.Button(self.menu_frame, text = "Double Peak",
                                    command=self.double_peak, state = tk.DISABLED)
        self.testbutton.pack(side=tk.LEFT)
        self.incbutton = ttk.Button(self.menu_frame, text = "+", command=self.increase_psize, state=tk.DISABLED)
        self.incbutton.pack(side=tk.LEFT)
        self.decbutton = ttk.Button(self.menu_frame, text = "-", command=self.decrease_psize, state=tk.DISABLED)
        self.decbutton.pack(side=tk.LEFT)
        self.runbutton = ttk.Button(self.menu_frame, text = "Run", command=self.run, state=tk.DISABLED)
        self.runbutton.pack(side=tk.LEFT)
        self.pizzabutton = ttk.Button(self.menu_frame, text = "Pizza!", command=self.makepizza, state=tk.NORMAL)
        self.pizzabutton.pack(side=tk.LEFT)
        self.configbutton = ttk.Button(self.menu_frame, text = "Configuration",
                                    command=self.config).pack(side=tk.LEFT)
        self.stylebutton = ttk.Button(self.menu_frame, text = "Theme",
                                    command = self.change_style).pack(side=tk.LEFT)
#Canvas event handling and other useful variables
        self.regex = r"\d{2,3}(?=\.)"
        self.peak_trigger = False
        self.run_trigger = False
        self.pizza_trigger = False
        self.peak_size = 6
        self.degree = 0
        self.colormap = colormaps.afmcold
        self.toggleViewState = tk.IntVar()
        self.toggleViewState.set(1)
        self.toggleShirleyState = tk.IntVar()
        self.toggleShirleyState.set(1)
        self.toggleGaussianState = tk.IntVar()
        self.toggleGaussianState.set(1)
        self.simmetryState = tk.IntVar()
        self.simmetryState.set(0)
#Define options for opening files
        self.file_opt = options = {}
        options['title'] = "Load your files"
        options['multiple'] = 1
        options['parent'] = self
        self.file_opt2 = options2 = {}
        self.files2 = []
        options2['title'] = "Load your output files"
        options2['parent'] = self
        options2['filetypes'] = [("easyPizza Output Files",("*.xpd","*.hol","*.dat"))]
        options2['initialdir'] = os.getcwd()
#State of the slide bar
        self.Scale = False
#Define the functions called
    def openfile(self):
        """Create a subplot and set index variable to 0 (first file).
        Open multiple text files with the data collected on LNLS, already splited.
        Then, activate the buttons that were Disabled and configure the Slider length."""
        self.ax = self.main_axes()
        self.index = 0
        self.files =  filedialog.askopenfilenames(**self.file_opt)
        if not self.files:
            return
        else:
            self.plotsimple(self.files[0])
            self.Scale = True
            self.evlabel["text"] = os.path.basename(self.files[0])
            self.slider.configure(from_=0,to=len(self.files)-1)
            self.testbutton["state"] = tk.NORMAL
            self.incbutton["state"] = tk.NORMAL
            self.decbutton["state"] = tk.NORMAL
            self.runbutton["state"] = tk.NORMAL
            if len(self.files) > 1:
                self.nextbutton["state"] = tk.NORMAL
                self.prevbutton["state"] = tk.NORMAL

    def plotsimple(self,txt):
        """The axis are cleared, in case any plot was there before.
        The regex searches for the energy in the text file, although this may be unnecessary.
        First the data is smoothed using a Savitzky Golay Filter, to prevent the Shirley to go nuts.
        There are options to use the Shirley Background Subtraction method, and adjust a Gaussian Peak.
        If data has two peaks, there is a Double Peak and + - buttons to adjust the size of the two peaks
        If the Run button is pressed, it integrates the Gaussian Peak adjusted and saves the data on output file."""
        self.ax.clear()
        self.ax.get_xaxis().set_ticks([])
        with open(txt) as File:
            L = re.findall(r"[-+]?\d*\.\d+|d+",File.readline())
        initialEnergy = float(L[0])
        finalEnergy = float(L[1])
        step = float(L[2])
        X = []
        for i in np.arange(initialEnergy,finalEnergy-0.1,step):
            X.append(i)

        Y = np.loadtxt(txt,skiprows=2)
        Y = savitzky.savitzky_golay(Y,7,3)

        if self.peak_trigger == False:
            if self.toggleShirleyState.get()==1:
                init, end = 0, len(X)

                Ys = shirley(Y,init,end-1)
                y_maxi = Ys.argmax()

                Xs = X-X[y_maxi]
                X = Xs
                Y = Y-min(Y)
                Xs, Ys = Xs[init+1:end], Ys[init+1:end]
                if self.toggleViewState.get() == 1:
                    self.ax.plot(Xs,Ys,'g-')
                if self.toggleGaussianState.get()==1:
                    p = gaussian.fit_iterative(Xs,Ys)
                    A, mu, sigma = p
                    x_n = np.linspace(np.min(Xs),np.max(Xs),100)
                    y = gaussian.gauss(x_n, A=A, mu=mu, sigma=sigma)
                    if self.toggleViewState.get() == 1:
                        self.ax.plot(x_n,y,'r-')
                    if self.run_trigger == True:
                        integrate = np.trapezoid(y,x_n)
                        txt = txt.split("/")[-1]
                        tp = re.findall(r"[-+]?\d*\.\d+|d+",txt)
                        tp = tp[0].split('.')
                        theta, phi = float(tp[0]), float(tp[1])
                        self.output.write("%f    %f    %f\n"%(integrate,theta,phi))
        else:
            if self.toggleShirleyState.get()==1:
                init, end = 0, len(X)

                Ys = shirley(Y,init,end-1)
                y_maxi = Ys.argmax()

                Xs = X-X[y_maxi]
                X = Xs
                Y = Y-min(Y)
                Xs, Ys = Xs[init+1:end], Ys[init+1:end]

                peak = max(Ys)
                peak_index = int(np.where(Ys==peak)[0][0])
                if peak_index>len(Xs)//2:
                    X1, Y1 = Xs[0:peak_index-self.peak_size], Ys[0:peak_index-self.peak_size]
                    X2, Y2 = Xs[peak_index-self.peak_size:-1], Ys[peak_index-self.peak_size:-1]
                else:
                    X1, Y1 = Xs[0:peak_index+self.peak_size], Ys[0:peak_index+self.peak_size]
                    X2, Y2 = Xs[peak_index+self.peak_size:-1], Ys[peak_index+self.peak_size:-1]

                XX = [X1, X2]
                YY = [Y1, Y2]

                if self.toggleViewState.get() == 1:
                    self.ax.plot(Xs,Ys,'g-')
                if self.toggleGaussianState.get()==1:
                    for a in range(len(XX)):
                        p = gaussian.fit_iterative(XX[a],YY[a])
                        A, mu, sigma = p
                        x_n = np.linspace(np.min(Xs),np.max(Xs),100)
                        y = gaussian.gauss(x_n, A=A, mu=mu, sigma=sigma)
                        if self.toggleViewState.get() == 1:
                            self.ax.plot(x_n,y,'r-')
                        if self.run_trigger == True:
                            integrate = np.trapezoid(y,x_n)
                            txt = txt.split("/")[-1]
                            tp = re.findall(r"[-+]?\d*\.\d+|d+",txt)
                            tp = tp[0].split('.')
                            theta, phi = float(tp[0]), float(tp[1])
                            if a==0:
                                self.output1.write("%f    %f    %f\n"%(integrate,theta,phi))
                            if a==1:
                                self.output2.write("%f    %f    %f\n"%(integrate,theta,phi))
        if self.toggleViewState.get() == 1:
            self.ax.plot(X,Y,'b-')
            self.ax.set_ylim(bottom=0)
            self.canvas.draw()

    def double_peak(self):
        """Activates or deactivates Double Peak"""
        if self.peak_trigger == False:
            self.peak_trigger = True
        else:
            self.peak_trigger = False
        self.plotsimple(self.files[self.index])

    def nextimg(self):
        """Go to the next plot"""
        if self.index < len(self.files)-1:
            self.index += 1
            self.slider.set(self.index)
            self.plotsimple(self.files[self.index])
            self.evlabel["text"] = os.path.basename(self.files[self.index])
            self.label["text"] = ""
        else:
            self.label["text"] = "Reached last image"

    def previmg(self):
        """Go to the previous plot"""
        if self.index == 0:
            self.label["text"] = "Reached first image"
        else:
            self.index -= 1
            self.slider.set(self.index)
            self.plotsimple(self.files[self.index])
            self.evlabel["text"] = os.path.basename(self.files[self.index])
            self.label["text"] = ""

    def increase_psize(self):
        """Increases the peak size in case of double peak."""
        self.peak_size += 1
        self.plotsimple(self.files[self.index])

    def decrease_psize(self):
        """Decreases the peak size in case of double peak."""
        self.peak_size -= 1
        self.plotsimple(self.files[self.index])
 
    def slide_image(self,value):
        """Uses the slider bar to move between plots"""
        if self.Scale:
            self.index = int(value)
            self.plotsimple(self.files[self.index])
            self.evlabel["text"] = os.path.basename(self.files[self.index])
            self.label["text"] = ""

    def config(self):
        """Sets a configuration window."""
        self.config_window = tk.Toplevel(self)
        self.config_frame = tk.Frame(self.config_window)
        self.config_frame.pack(side='top',expand=True, fill='both')
        self.barfig = Figure(figsize=(5,3))
        self.barfig.set_facecolor("#d9d9d9")
        self.ax1 = self.barfig.add_axes([0.1,0,0.1,0.9])
        self.ax2 = self.barfig.add_axes([0.3,0,0.1,0.9])
        self.ax3 = self.barfig.add_axes([0.6,0,0.1,0.9])
        self.ax4 = self.barfig.add_axes([0.8,0,0.1,0.9])
        cb1 = mpl.colorbar.ColorbarBase(self.ax1, cmap = colormaps.afmcold)
        cb2 = mpl.colorbar.ColorbarBase(self.ax2, cmap = colormaps.afmgreen)
        cb3 = mpl.colorbar.ColorbarBase(self.ax3, cmap = colormaps.afmhot)
        cb4 = mpl.colorbar.ColorbarBase(self.ax4, cmap = colormaps.inferno)
        cb1.ax.set_yticklabels([])
        cb2.ax.set_yticklabels([])
        cb3.ax.set_yticklabels([])
        cb4.ax.set_yticklabels([])
        self.barcanvas = FigureCanvasTkAgg(self.barfig,master=self.config_frame)
        self.barcanvas.get_tk_widget().pack(side='top', fill=tk.BOTH, expand=True)
        self.barcanvas.draw()
        self.button1 = tk.Button(self.config_frame, text="Blue",command=self.select_blue,width=10)
        self.button2 = tk.Button(self.config_frame, text="Green",command=self.select_green,width=10)
        self.button3 = tk.Button(self.config_frame, text="Red",command=self.select_red,width=10)
        self.button4 = tk.Button(self.config_frame, text="Inferno",command=self.select_inferno,width=10)
        self.button1.pack(side='left')
        self.button2.pack(side='left')
        self.button3.pack(side='left')
        self.button4.pack(side='left')
        self.barcanvas.draw()
        self.checkV_frame = tk.Frame(self.config_window)
        self.checkV_frame.pack(side='bottom',expand=True,fill='both')
        self.checkS_frame = tk.Frame(self.config_window)
        self.checkS_frame.pack(side='bottom',expand=True,fill='both')
        self.checkG_frame = tk.Frame(self.config_window)
        self.checkG_frame.pack(side='bottom',expand=True,fill='both')
        self.toggleView = tk.Checkbutton(self.checkV_frame,variable=self.toggleViewState, onvalue=1,
                                            offvalue=0, text="Toggle graph visualization").pack(side='left')
        self.toggleShirley = tk.Checkbutton(self.checkS_frame,variable=self.toggleShirleyState, onvalue=1,
                                            offvalue=0, text ="Shirley Background Subtraction").pack(side='left')
        self.toggleGaussian = tk.Checkbutton(self.checkG_frame,variable=self.toggleGaussianState, onvalue=1,
                                            offvalue=0, text ="Gaussian Fit").pack(side='left')
        self.pizzaframe = tk.Frame(self.config_window)
        self.pizzaframe.pack(side="top", expand=True, fill="both")
        self.rotatelabel = tk.Label(self.pizzaframe, text="Rotate Pizza (degrees)").pack()
        self.r1 = tk.Entry(self.pizzaframe)
        self.r1.pack(padx=5)
        self.rotatebutton = tk.Button(self.pizzaframe, text="Rotate", command=self.rotate)
        self.rotatebutton.pack()
        self.simmetryframe = tk.Frame(self.config_window)
        self.simmetryframe.pack(side="top", expand=True, fill="both")
        self.simmetrylabel = tk.Label(self.simmetryframe, text="Select system simmetry").pack()
        self.raw = tk.Radiobutton(self.simmetryframe, text = "Raw Data",
                    variable=self.simmetryState, value=0,command=self.change_simmetry).pack()
        self.d60 = tk.Radiobutton(self.simmetryframe, text = u"60\N{DEGREE SIGN}", 
                    variable = self.simmetryState, value=1, command = self.change_simmetry).pack()
        self.d120 = tk.Radiobutton(self.simmetryframe, text = u"120\N{DEGREE SIGN}", 
                    variable = self.simmetryState, value=2, command = self.change_simmetry).pack()
        self.d90 = tk.Radiobutton(self.simmetryframe, text = u"90\N{DEGREE SIGN}",
                    variable = self.simmetryState, value=3, command = self.change_simmetry).pack()
        self.exportframe = tk.Frame(self.config_window)
        self.exportframe.pack(side="top", expand=True, fill="both")
        self.exportbutton = tk.Button(self.exportframe, text="Export Data", command=self.export)
        self.exportbutton.pack()

    def export(self):
        try:
            a = self.files2.split("/")[-1].split(".")
            for tt in range(int(min(self.THETA)),int(max(self.THETA))+1,3):
                for pp in range(int(360*min(self.PHI)/(2*np.pi)),int(360*max(self.PHI)/(2*np.pi))+1,3):
                    p = self.find_nearest(self.gphi*360/(2*np.pi),pp)
                    t = self.find_nearest(self.gtheta,tt)
                    print(pp, tt)
        except NameError:
            pass

    def find_nearest(self, array, value):
        idx = (np.abs(array-value)).argmin()
        idx2 = np.unravel_index(idx,array.shape)
        return idx2

    def change_simmetry(self):
        if self.simmetryState.get()==0:
            self.plotpizza(rotate=self.degree)
        if self.simmetryState.get()==1:
            self.plotpizza(rotate=self.degree)
        if self.simmetryState.get()==2:
            self.plotpizza(rotate=self.degree)
        if self.simmetryState.get()==3:
            self.plotpizza(rotate=self.degree)

    def select_blue(self,event=None):
        """Selects the afmcold colormap."""
        self.colormap = colormaps.afmcold
        self.plotpizza(rotate = self.degree)

    def select_green(self,event=None):
        """Selects the afmgreen colormap."""
        self.colormap = colormaps.afmgreen
        self.plotpizza(rotate = self.degree)

    def select_red(self,event=None):
        """Selects the afmhot colormap."""
        self.colormap = colormaps.afmhot
        self.plotpizza(rotate = self.degree)

    def select_inferno(self,event=None):
        """Selects the inferno colormap."""
        self.colormap = colormaps.inferno
        self.plotpizza(rotate = self.degree)

    def change_style(self):
        """Set a new window with a Combobox to choose between themes"""
        self.combo_window = tk.Toplevel(self)
        self.combo = ttk.Combobox(self.combo_window, values=self.available_themes)
        self.combo.pack(padx=32, pady=8)
        self.combo.bind('<Return>', self.set_style)
        button = ttk.Button(self.combo_window, text='OK')
        button['command'] = self.set_style
        button.pack(pady=8)

    def set_style(self, event=None):
        """Set the Style to the content of the Combobox"""
        content = self.combo.get()
        try:
            self.style.theme_use(content)
            self.style.configure('TButton', relief='flat')
        except tk.TclError as err:
            messagebox.showerror('Error', err)
        self.combo_window.destroy()

    def run(self):
        """Asks for the name of the output file(s) and go through every plot, 
        integrating the data and saving to the outputs."""
        self.toggleShirleyState.set(1)
        self.toggleGaussianState.set(1)
        self.out_window = tk.Toplevel(self)
        if self.peak_trigger == True:
            tk.Label(self.out_window, text="Output File 1").pack()
            self.e1 = tk.Entry(self.out_window)
            self.e1.pack(padx=5)
            tk.Label(self.out_window, text="Output File2").pack()
            self.e2 = tk.Entry(self.out_window)
            self.e2.pack(padx=5)
        else:
            tk.Label(self.out_window, text="Output File").pack()
            self.e = tk.Entry(self.out_window)
            self.e.pack(padx=5)
        self.okbutton = tk.Button(self.out_window, text="OK",command=self.ok)
        self.okbutton.pack(pady=5)
        self.out_window.bind('<Return>', self.ok)

    def ok(self,event=None):
        """Gets the state of the Label on the Run window and create the output files.
        Creates progress bar to see remaining time."""
        self.run_trigger = True
        self.index = 0
        if self.peak_trigger == True:
            self.output1 = open(self.e1.get()+".xpd",'w')
            self.output1.write("Intensity    Theta    Phi\n")
            self.output2 = open(self.e2.get()+".xpd",'w')
            self.output2.write("Intensity    Theta    Phi\n")
        else:
            self.output = open(self.e.get()+".xpd",'w')
            self.output.write("Intensity    Theta    Phi\n")
        self.out_window.destroy()
        self.progress["value"] = 0
        self.progress["maximum"] = len(self.files)-1
        for i in range(len(self.files)):
            self.plotsimple(self.files[self.index])
            self.evlabel["text"] = os.path.basename(self.files[self.index])
            self.index += 1
            self.progress["value"] = i
            self.update_idletasks()
        if self.peak_trigger==True:
            self.output1.close()
            self.output2.close()
        else:
            self.output.close()
        self.run_trigger = False

    def rotate(self):
        """Rotates Pizza"""
        a = self.r1.get()
        try:
            self.degree = float(a)
            self.plotpizza(rotate=self.degree)
        except ValueError:
            pass


    def makepizza(self):
        """Asks for an output file (searches for a .xpd file, generated in the program.
        The .xpd file must be a file with Intensity, Theta, Phi, in this order.
        Generates the XPD pizza, in this case with 120 degree simmetry, although it can be changed easilly."""
        if self.pizza_trigger==True:
            self.config_window.destroy()
            self.pizza_trigger = False
        self.config()
        self.files2 =  filedialog.askopenfilename(**self.file_opt2)
        self.plotpizza()
        self.pizza_trigger = True

    def main_axes(self):
        """Return the Cartesian axes, creating it only if it is no longer in the figure.
        Old Matplotlib returned the existing axes from add_subplot(111); new versions stack a new one."""
        if self.ax not in self.fig.axes:
            self.ax = self.fig.add_subplot(111)
        return self.ax

    def plotpizza(self, rotate=0):
        self.ax = self.main_axes()
        self.ax.clear()
        if not self.files2:
            return
        symmetry = self.simmetryState.get()
        THETA, PHI, ki = angular_map.load_map(self.files2, symmetry, rotate)
        self.THETA, self.PHI = THETA, PHI
        self.gtheta, self.gphi, self.z = angular_map.interpolate_map(THETA, PHI, ki, symmetry)

        norm = mpl.colors.Normalize(self.z.flatten().min(),self.z.flatten().max())

        self.fig.delaxes(self.ax)
        self.ax2 = self.fig.add_subplot(1,1,1, polar = True)
        # Matplotlib < 3.3 drew a same-shape C with 'flat' shading by dropping its last row and column.
        self.ax2.pcolormesh(self.gphi, self.gtheta, self.z[:-1, :-1], shading='flat', cmap=self.colormap, norm=norm)
        self.ax2.set_yticks([min(THETA),max(THETA)])
        self.ax2.grid(True)
        self.ax2.set_rmax(max(THETA))
        self.fig.savefig(self.files2+".jpg")
        self.canvas.draw()
        self.fig.delaxes(self.ax2)

if __name__ == "__main__":
    MainWindow().mainloop()
