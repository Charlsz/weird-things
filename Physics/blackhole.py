#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Black Hole Photon Orbit Visualizer with PyQt5 Interface
Real-time visualization of null geodesics in Schwarzschild spacetime
with accretion disk rendering
"""

import warnings
warnings.filterwarnings('ignore')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys
import matplotlib
matplotlib.use('Qt5Agg')
from PyQt5 import QtCore, QtWidgets, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from scipy.integrate import odeint
import seaborn as sns
from scipy.interpolate import interp1d
from scipy.special import ellipj, ellipkinc, ellipk
from scipy.optimize import fsolve

# ============================================================
#  NULL GEODESIC EQUATIONS (Light Ray Tracing)
# ============================================================

def S_null(Z, t, p):
    """
    System of ODEs for null geodesics in Schwarzschild metric
    Z = [r, rdot, phi]
    p = [M, L] where M is mass and L is angular momentum
    """
    r, rdot, phi = Z
    M, L = p
    phidot = L / r**2
    return [
        rdot,
        L**2 * (r - 3*M) / r**4,
        phidot
    ]


def init_cond(b, x_init):
    """
    Initial conditions for a photon with impact parameter b
    starting at position x_init
    """
    r_init = np.sqrt(b**2 + x_init**2)
    phi_init = np.arccos(x_init / r_init)
    rdot_init = np.cos(phi_init)
    phidot_init = -np.sqrt((1 - rdot_init**2) / r_init**2)
    L = r_init**2 * phidot_init
    return [r_init, rdot_init, phi_init], L


def integrate(t, initial, p):
    """
    Integrate the geodesic equations
    """
    sol = odeint(S_null, initial, t, args=(p,))
    r = sol[:, 0]
    phi = sol[:, 2]
    return r, phi


# ============================================================
#  ACCRETION DISK ISORADIAL CURVES (Elliptic Functions)
# ============================================================

def B_fun(p, M):
    """Impact parameter as function of orbital radius"""
    return (p**3 / (p - 2*M))**0.5


def Q_fun(P, M):
    """Q parameter for elliptic integrals"""
    return ((P - 2*M) * (P + 6*M))**0.5


def gamma(alpha, theta_0):
    """Angle transformation for observer inclination"""
    a = np.cos(alpha)**2 + 1 / (np.tan(theta_0)**(2))
    return np.arccos(np.cos(alpha) / (a**0.5))


def k2(P, M):
    """Elliptic modulus squared"""
    Q = Q_fun(P, M)
    return ((Q - P + 6*M) / (2*Q))


def zeta_inf(P, M):
    """Asymptotic angle for elliptic integral"""
    Q = Q_fun(P, M)
    ratio = (Q - P + 2*M) / (Q - P + 6*M)
    return np.arcsin(ratio**0.5)


def Up(P, alpha, M, theta_0):
    """
    First order image: 1/r as function of angle
    Using Jacobi elliptic functions
    """
    Q = Q_fun(P, M)
    A1 = (Q - P + 2*M) / (4*M*P)
    A2 = (Q - P + 6*M) / (4*M*P)
    g = gamma(alpha, theta_0)
    ratio = (Q / P)**0.5
    mod = k2(P, M)
    sn, cn, dn, ph = ellipj(((g/2)*ratio + ellipkinc(zeta_inf(P, M), mod)), mod)
    return -A1 + A2*sn**2


def Up2(P, alpha, M, theta_0):
    """
    Second order image: photons that circle the BH before reaching observer
    """
    Q = Q_fun(P, M)
    A1 = (Q - P + 2*M) / (4*M*P)
    A2 = (Q - P + 6*M) / (4*M*P)
    g = gamma(alpha, theta_0)
    ratio = (Q / P)**0.5
    mod = k2(P, M)
    sn, cn, dn, ph = ellipj((0.5*(g-2*np.pi)*ratio + 2*ellipk(mod) - 
                             ellipkinc(zeta_inf(P, M), mod)), mod)
    return -A1 + A2*sn**2


# ============================================================
#  PLOTTING FUNCTION
# ============================================================

def plot_solution(self, fig, x_1, y_1, xs_acc_disk, ys_acc_disk, 
                 xs_acc_disk_2nd, ys_acc_disk_2nd, angle_disk, 
                 accret_disk_min_r, accret_disk_max_r, M, t):
    """
    Create the two-panel visualization:
    Left: Geodesics with accretion disk
    Right: Observer view of accretion disk
    """
    
    # Left panel: Geodesics
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.grid(False)
    
    # Draw accretion disk line
    x_disk_plus = np.arange(accret_disk_min_r, accret_disk_max_r, 0.05)
    x_disk_minus = np.arange(-accret_disk_max_r, -accret_disk_min_r, 0.05)
    
    # Check which photons hit the disk
    for i in range(len(x_1)):
        attempts = 0
        hit_disk = False
        
        # Check positive side
        for j in range(len(x_disk_plus)):
            attempts += 1
            d = np.sqrt((x_1[i] - x_disk_plus[j])**2 + 
                       (y_1[i] - x_disk_plus[j]*np.tan(angle_disk))**2)
            if np.min(d) < 0.1:
                idx = j
                plt.scatter(x_disk_plus[idx], x_disk_plus[idx]*np.tan(angle_disk), 
                          color='black', s=20)
                line_hit, = plt.plot(x_1[i], y_1[i], color='blue', linewidth=1.5)
                hit_disk = True
                break
        
        # Check negative side
        attempts_2 = 0
        if not hit_disk:
            for j in range(len(x_disk_minus)):
                attempts_2 += 1
                d = np.sqrt((x_1[i] - x_disk_minus[j])**2 + 
                           (y_1[i] - x_disk_minus[j]*np.tan(angle_disk))**2)
                if np.min(d) < 0.1:
                    idx = j
                    plt.scatter(x_disk_minus[idx], x_disk_minus[idx]*np.tan(angle_disk), 
                              color='black', s=20)
                    plt.plot(x_1[i], y_1[i], color='blue', linewidth=1.5)
                    hit_disk = True
                    break
        
        # If didn't hit disk
        if not hit_disk:
            line_fail, = plt.plot(x_1[i], y_1[i], color='lightgrey', linewidth=1)
    
    # Draw black hole horizon
    circle = plt.Circle((0., 0.), 2*M, color='black', fill=True, zorder=10)
    plt.gca().add_patch(circle)
    
    # Draw accretion disk
    plt.plot(x_disk_plus, x_disk_plus*np.tan(angle_disk), color='red', linewidth=2)
    plt.plot(x_disk_minus, x_disk_minus*np.tan(angle_disk), color='red', linewidth=2)
    
    # Labels and styling
    plt.xlabel(r"Distance $X$", fontsize=10, color='white')
    plt.ylabel(r"Distance $Z$", fontsize=10, color='white')
    plt.yticks(color='white')
    plt.xticks(color='white')
    plt.axis('equal')
    
    # Legend
    try:
        ax1.legend(handles=[line_hit, line_fail], 
                  labels=['Accretion disk hit', 'Off to infinity or captured'], 
                  loc='upper right')
    except UnboundLocalError:
        print('All lightrays either hit the accretion disk or miss it.')
    
    ax1.set_xlim(-20, 20)
    ax1.set_xbound(-20, 20)
    ax1.set_ylim(-20, 20)
    ax1.set_ybound(-20, 20)
    ax1.set_title(r"Lightlike geodesics for" "\n" 
                 r"$ds^2 = -(1-2M/r) \ dt^2 + (1 - 2M/r)^{-1}dr^2 + r^2 d\Omega^2$",
                 fontsize=12, color='white')
    
    # Right panel: Accretion disk view
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.grid(False)
    
    # Draw isoradial curves with hot colormap
    for i in range(len(xs_acc_disk)):
        x = xs_acc_disk[i]
        y = ys_acc_disk[i]
        x2 = xs_acc_disk_2nd[i]
        y2 = ys_acc_disk_2nd[i]
        
        # Interpolate for smooth curves at high inclinations
        if np.pi/2 - angle_disk < 78*np.pi/180:
            if i > 3 and len(y) > 40:
                try:
                    x_interp = np.concatenate((y[-20:-10], y[10:20]))
                    y_interp = np.concatenate((-x[-20:-10], -x[10:20]))
                    f = interp1d(x_interp, y_interp, kind='quadratic')
                    x_new = np.linspace(y[-10], y[10], 100)
                    plt.plot(x_new, f(x_new), 
                           color=sns.color_palette('hot_r', len(xs_acc_disk))[i],
                           zorder=2, lw=4)
                except:
                    pass
        
        # Second order images
        if len(y2) > 0 and len(x2) > 0:
            plt.plot(y2, x2, color=sns.color_palette('hot_r', len(xs_acc_disk))[i],
                    zorder=1, lw=4)
        
        # First order images
        if len(y) > 0 and len(x) > 0:
            plt.plot(y, -x, color=sns.color_palette('hot_r', len(xs_acc_disk))[i],
                    zorder=2, lw=4)
    
    # Black hole and photon sphere
    circle = plt.Circle((0., 0.), 2*M, color='black', fill=True, zorder=100)
    plt.gca().add_patch(circle)
    circle = plt.Circle((0., 0.), 3*M, color='black', fill=True, zorder=99, alpha=0.3)
    plt.gca().add_patch(circle)
    circle = plt.Circle((0., 0.), np.sqrt(3)*3*M, facecolor='none', 
                       edgecolor='white', fill=False, zorder=98, lw=4)
    plt.gca().add_patch(circle)
    
    ax2.set_xlim(-32, 32)
    ax2.set_ylim(-32, 32)
    ax2.axis('off')
    fig.patch.set_color('dimgrey')
    for ax in fig.axes:
        ax.patch.set_color('dimgrey')
    ax2.set_title(f'Black hole accretion disk view at {int((np.pi/2 - angle_disk)*180/np.pi)}º', 
                 fontsize=14, color='white')
    
    print('Done! Check your plots :)')


# ============================================================
#  PyQt5 GUI WIDGET
# ============================================================

class Widget(QtWidgets.QWidget):
    def __init__(self, *args, **kwargs):
        QtWidgets.QWidget.__init__(self, *args, **kwargs)
        
        # Create matplotlib figure
        global fig
        fig = plt.figure(figsize=(12, 6), dpi=100)
        fig.set_facecolor('#1a1a1a')
        self.figure = fig
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QtWidgets.QSizePolicy.Expanding, 
                                 QtWidgets.QSizePolicy.Expanding)
        
        # Main layout
        main_layout = QtWidgets.QVBoxLayout(self)
        
        # Title
        title = QtWidgets.QLabel("BLACK HOLE SIMULATOR")
        title.setStyleSheet("""
            color: white;
            font-size: 24px;
            font-weight: bold;
            padding: 15px;
            background-color: #2a2a2a;
            border-radius: 5px;
        """)
        title.setAlignment(QtCore.Qt.AlignCenter)
        main_layout.addWidget(title)
        
        # Canvas
        main_layout.addWidget(self.canvas)
        
        # Simple controls
        control_panel = QtWidgets.QWidget()
        control_panel.setStyleSheet("background-color: #2a2a2a; padding: 10px; border-radius: 5px;")
        control_layout = QtWidgets.QGridLayout()
        
        # Number of photons
        photon_label = QtWidgets.QLabel("Number of Photons:")
        photon_label.setStyleSheet("color: white; font-size: 14px;")
        self.photon_spin = QtWidgets.QSpinBox()
        self.photon_spin.setRange(1, 50)
        self.photon_spin.setValue(10)
        self.photon_spin.setStyleSheet("""
            QSpinBox {
                background-color: #3a3a3a;
                color: white;
                font-size: 14px;
                padding: 5px;
                border: 1px solid #555;
                border-radius: 3px;
            }
        """)
        
        control_layout.addWidget(photon_label, 0, 0)
        control_layout.addWidget(self.photon_spin, 0, 1)
        
        # View angle
        angle_label = QtWidgets.QLabel("View Angle (degrees):")
        angle_label.setStyleSheet("color: white; font-size: 14px;")
        self.angle_spin = QtWidgets.QSpinBox()
        self.angle_spin.setRange(0, 85)
        self.angle_spin.setValue(20)
        self.angle_spin.setStyleSheet("""
            QSpinBox {
                background-color: #3a3a3a;
                color: white;
                font-size: 14px;
                padding: 5px;
                border: 1px solid #555;
                border-radius: 3px;
            }
        """)
        
        control_layout.addWidget(angle_label, 1, 0)
        control_layout.addWidget(self.angle_spin, 1, 1)
        
        # Buttons
        btn_run = QtWidgets.QPushButton("RUN SIMULATION")
        btn_run.clicked.connect(self.main)
        btn_run.setStyleSheet("""
            QPushButton {
                background-color: #4a4a4a;
                color: white;
                font-size: 14px;
                font-weight: bold;
                padding: 10px 20px;
                border: 1px solid #666;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #5a5a5a;
            }
        """)
        
        btn_stop = QtWidgets.QPushButton("STOP")
        btn_stop.clicked.connect(self.stop)
        btn_stop.setStyleSheet("""
            QPushButton {
                background-color: #4a4a4a;
                color: white;
                font-size: 14px;
                font-weight: bold;
                padding: 10px 20px;
                border: 1px solid #666;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #5a5a5a;
            }
        """)
        
        control_layout.addWidget(btn_run, 2, 0)
        control_layout.addWidget(btn_stop, 2, 1)
        
        control_panel.setLayout(control_layout)
        main_layout.addWidget(control_panel)
        
        self.setStyleSheet("background-color: #1a1a1a;")
    
    def case(self, text):
        """Handle predefined cases"""
        if text == "Critical Loop":
            self.Loop()
    
    def stop(self):
        """Stop the current simulation"""
        try:
            global ani
            ani.event_source.stop()
            ani.frame_seq = ani.new_frame_seq()
        except:
            pass
    
    def Loop(self):
        """
        Predefined case: critical impact parameter
        Shows photon orbiting many times before capture/escape
        """
        self.figure.clear()
        M = 1
        angle_disk = 0
        theta_0 = np.pi/2
        accret_disk_min_r = 6*M
        accret_disk_max_r = 15*M
        r_list = np.arange(4, 30, 0.5)*M
        b_c = 3*np.sqrt(3)*M
        xs_acc_disk = []
        ys_acc_disk = []
        xs_acc_disk_2nd = []
        ys_acc_disk_2nd = []
        
        print("Computing isoradial curves... (this may take a moment)")
        
        # Compute isoradial curves (simplified for speed)
        for idx, r in enumerate(r_list[::2]):  # Skip some for speed
            alpha_list = np.arange(0, 4*np.pi, 0.05)
            b_list = []
            alpha_res = []
            
            for alpha in alpha_list:
                def cu(P, M, theta_0):
                    return 1 - r*Up(P, alpha, M, theta_0)
                
                guess = b_c + 0.1
                try:
                    raiz = fsolve(cu, [guess], args=(M, theta_0), full_output=False)[0]
                    if abs(raiz - guess) > 0.01:
                        b_raiz = B_fun(raiz, M)
                        if b_raiz > 0 and b_raiz < 100:
                            b_list.append(b_raiz)
                            alpha_res.append(alpha)
                except:
                    pass
            
            if len(b_list) > 0:
                b_array = np.array(b_list)
                alpha_array = np.array(alpha_res)
                x = b_array * np.cos(alpha_array)
                y = b_array * np.sin(alpha_array)
                xs_acc_disk.append(x)
                ys_acc_disk.append(y)
                xs_acc_disk_2nd.append(np.array([]))
                ys_acc_disk_2nd.append(np.array([]))
        
        # Light ray at critical impact parameter
        x_init = -40
        max_t = 100.
        b_min = 3*np.sqrt(3)*M
        b_max = 3*np.sqrt(3)*M
        b_num = 1
        b_list = np.linspace(b_min, b_max, int(b_num))
        t = np.arange(0, max_t, 0.01)
        
        xs = []
        ys = []
        for j, b in enumerate(b_list):
            initial, L = init_cond(b, x_init)
            r, phi = integrate(t, initial, [M, L])
            
            x_1 = r * np.cos(phi)
            y_1 = r * np.sin(phi)
            xs.append(x_1)
            ys.append(y_1)
        
        plot_solution(self, fig, xs, ys, xs_acc_disk, ys_acc_disk, 
                     xs_acc_disk_2nd, ys_acc_disk_2nd, angle_disk, 
                     accret_disk_min_r, accret_disk_max_r, M, t)
        self.canvas.draw()
    
    def main(self):
        """
        Main simulation with user-defined parameters
        """
        self.figure.clear()
        M = 1
        angle_disk = self.angle_spin.value() * np.pi / 180
        theta_0 = np.pi/2 - self.angle_spin.value() * np.pi / 180
        accret_disk_min_r = 6*M
        accret_disk_max_r = 15*M
        r_list = np.arange(4, 30, 0.5)*M
        b_c = 3*np.sqrt(3)*M
        xs_acc_disk = []
        ys_acc_disk = []
        xs_acc_disk_2nd = []
        ys_acc_disk_2nd = []
        
        print("Computing accretion disk isoradial curves...")
        
        # Compute isoradial curves
        for idx, r in enumerate(r_list[::2]):  # Skip for performance
            if theta_0 < 78*np.pi/180:
                alpha_list = np.arange(0, 2*np.pi, 0.05)
            else:
                alpha_list = np.arange(0, 4*np.pi, 0.05)
            
            b_list = []
            b_list2 = []
            alpha_res = []
            alpha_res2 = []
            
            for alpha in alpha_list:
                # First order
                def cu(P, M, theta_0):
                    return 1 - r*Up(P, alpha, M, theta_0)
                
                # Second order
                def cu2(P, M, theta_0):
                    return -(1 - r*Up2(P, alpha, M, theta_0))
                
                guess = b_c + 0.1
                
                try:
                    raiz = fsolve(cu, [guess], args=(M, theta_0))[0]
                    if abs(raiz - guess) > 0.01:
                        b_raiz = B_fun(raiz, M)
                        if b_raiz > 0 and b_raiz < 100:
                            b_list.append(b_raiz)
                            alpha_res.append(alpha)
                except:
                    pass
                
                try:
                    raiz2 = fsolve(cu2, [guess], args=(M, theta_0))[0]
                    if abs(raiz2 - guess) > 0.01:
                        b_raiz2 = B_fun(raiz2, M)
                        if b_raiz2 > 0 and b_raiz2 < 100:
                            b_list2.append(b_raiz2)
                            alpha_res2.append(alpha)
                except:
                    pass
            
            # First order
            if len(b_list) > 0:
                b_array = np.array(b_list)
                alpha_array = np.array(alpha_res)
                x = b_array * np.cos(alpha_array)
                y = b_array * np.sin(alpha_array)
                xs_acc_disk.append(x)
                ys_acc_disk.append(y)
            else:
                xs_acc_disk.append(np.array([]))
                ys_acc_disk.append(np.array([]))
            
            # Second order
            if len(b_list2) > 0:
                b_array2 = np.array(b_list2)
                alpha_array2 = np.array(alpha_res2)
                x2 = b_array2 * np.cos(alpha_array2)
                y2 = b_array2 * np.sin(alpha_array2)
                
                # Filter points within view limits
                xlims = [-32, 32]
                ylims = [-32, 32]
                mask = (x2 > xlims[0]) & (x2 < xlims[1]) & (y2 > ylims[0]) & (y2 < ylims[1])
                x2_filtered = x2[mask]
                y2_filtered = y2[mask]
                
                xs_acc_disk_2nd.append(x2_filtered)
                ys_acc_disk_2nd.append(y2_filtered)
            else:
                xs_acc_disk_2nd.append(np.array([]))
                ys_acc_disk_2nd.append(np.array([]))
        
        print("Tracing light rays...")
        
        # Trace light rays
        x_init = -40
        max_t = 100.
        critical = 3*np.sqrt(3)*M
        b_min = critical - 2
        b_max = critical + 2
        b_num = self.photon_spin.value()
        b_list = np.linspace(b_min, b_max, int(b_num))
        t = np.arange(0, max_t, 0.01)
        
        xs = []
        ys = []
        for j, b in enumerate(b_list):
            initial, L = init_cond(b, x_init)
            r, phi = integrate(t, initial, [M, L])
            
            if b < 0:
                x_1 = r * np.cos(phi)
                y_1 = -r * np.sin(phi)
            else:
                x_1 = r * np.cos(phi)
                y_1 = r * np.sin(phi)
            xs.append(x_1)
            ys.append(y_1)
        
        plot_solution(self, fig, xs, ys, xs_acc_disk, ys_acc_disk, 
                     xs_acc_disk_2nd, ys_acc_disk_2nd, angle_disk, 
                     accret_disk_min_r, accret_disk_max_r, M, t)
        self.canvas.draw()


# ============================================================
#  MAIN APPLICATION
# ============================================================

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle('Fusion')
    
    w = Widget()
    w.setWindowTitle('Black Hole Simulator')
    w.resize(1400, 800)
    w.show()
    
    print("Black Hole Simulator")
    print("Adjust number of photons and viewing angle, then click RUN SIMULATION\n")
    
    sys.exit(app.exec_())
    w.show()
    sys.exit(app.exec_())
