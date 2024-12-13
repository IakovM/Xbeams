# For details of Crossed molecular beams experiments please see:
# [1] Gu, X., Guo, Y., Zhang, F., Mebel, A.M. and Kaiser, R.I., 2006. Reaction dynamics of carbon-bearing radicals in circumstellar envelopes of carbon stars. Faraday discussions, 133, pp.245-275. 
# DOI	https://doi.org/10.1039/B516457E
# [2] Kaiser, R.I., 2002. Experimental investigation on the formation of carbon-bearing molecules in the interstellar medium via neutral− neutral reactions. Chemical Reviews, 102(5), pp.1309-1358.
# https://doi.org/10.1021/cr970004v


def center(a_mass, b_mass, c_mass, d_mass, a_vp, a_s, b_vp, b_s, offset, er, A, B, C, D, A_carrier_gas, x_max=2200, x_min=-200, y_max=2200, y_min=-200, save_pdf=False):
    #This function solves the center-of-mass (CM) angle, Ec and so on (collisional energy, max velocity, relative velocity in CM).
    #a_mass, mass of a
     #b_mass,
     #c_mass,
     #d_mass,
     #a_vp, most prob vel of a
     #a_s, S of A beam
     #b_vp, most prob vel of B
     #b_s, S of B Beam
     #offset, offset to calculate TOF parameters
     #er, reaction energy
     #A, name of A
     #B, name of B
     #C, name of C
     #D, name of C
     #A_carrier_gas,
     #x_max = 2200, x on for CM fig
    #x_min = -200, x on for CM fig
    #y_max = 2200, y on for CM fig
    #y_min = -200 y on for CM fig
    # save_pdf=False = make True if you want to save a pdf with the name CM.pdf
    import math
    import matplotlib.pyplot as plt
    import numpy as np

    # For future LaTeX not to be ugly
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Arial'
    plt.rcParams['mathtext.it'] = 'Arial'
    plt.rcParams['mathtext.bf'] = 'Arial'
    plt.rc('font', family='arial', weight='normal', size=12)

    # offset of the detector
    cm_offset = 0.75

    # ec = energy of collision
    ec = 0.5 * a_mass * b_mass * 0.001 * (a_vp * a_vp + b_vp * b_vp) / (a_mass + b_mass)

    x = ((er + ec) * 2 * 1000 * (a_mass + b_mass) / c_mass / d_mass) ** 0.5
    # max recoil velocity of a heavy product aka n_radius
    n_radius = d_mass * x / (d_mass + c_mass)
    #speed in CM frame
    cm_speed = (a_mass ** 2 * a_vp ** 2 + b_mass ** 2 * b_vp ** 2) ** 0.5 / (a_mass + b_mass)

    print(n_radius)
    print(cm_speed)

    if n_radius / cm_speed >= 1.0:
        high_limit = math.atan(b_mass * b_vp / a_mass / a_vp) * 180 / math.pi + cm_offset + math.asin(1) * 180 / math.pi
        low_limit = math.atan(b_mass * b_vp / a_mass / a_vp) * 180 / math.pi + cm_offset - math.asin(1) * 180 / math.pi
    else:
        high_limit = math.atan(b_mass * b_vp / a_mass / a_vp) * 180 / math.pi + cm_offset + math.asin(n_radius / cm_speed) * 180 / math.pi
        low_limit = math.atan(b_mass * b_vp / a_mass / a_vp) * 180 / math.pi + cm_offset - math.asin(n_radius / cm_speed) * 180 / math.pi

    cm = math.atan(b_mass * b_vp / a_mass / a_vp) * 180 / math.pi + cm_offset
    # high_limit = math.atan(b_mass*b_vp/a_mass/a_vp)*180/math.pi+cm_offset+math.asin(n_radius/cm_speed)*180/math.pi
    # low_limit =  math.atan(b_mass*b_vp/a_mass/a_vp)*180/math.pi+cm_offset-math.asin(n_radius/cm_speed)*180/math.pi
    x_cm = a_vp / (math.tan(math.pi / 2 - math.atan(b_mass * b_vp / a_mass / a_vp)) + a_vp / b_vp)
    y_cm = x_cm * math.tan(math.pi / 2 - math.atan(b_mass * b_vp / a_mass / a_vp))
    # low_limit_x = a_vp/(math.tan(math.pi/2-math.atan(b_mass*b_vp/a_mass/a_vp)+math.asin(n_radius/cm_speed))+a_vp/b_vp)
    # low_limit_y = low_limit_x*math.tan(math.pi/2-math.atan(b_mass*b_vp/a_mass/a_vp)+math.asin(n_radius/cm_speed))
    # high_limit_x = a_vp
    # high_limit_y = high_limit_x*math.tan(math.pi/2-(math.atan(b_mass*b_vp/a_mass/a_vp)+math.asin(n_radius/cm_speed)))
    # high_limit_y

    # calculate v_cm of a product
    v_cm = math.sqrt(x_cm ** 2 + y_cm ** 2)

    # estimate time of flight of the product (parameters are specific for the UH Manoa XB machine)
    t_o_f = 329 / 1000 / v_cm * 10 ** 6 + offset
    t_o_f_ch = t_o_f / 10.24
    t_o_f_low = 329 / 1000 / (v_cm + n_radius) * 10 ** 6 + offset
    t_o_f_low_ch = t_o_f_low / 10.24
    t_o_f_high = 329 / 1000 / (v_cm - n_radius) * 10 ** 6 + offset
    t_o_f_high_ch = t_o_f_high / 10.24

    # now plot a Single reactive channel

    plt.figure(figsize=(10, 10))
    plt.plot([0, 0], [0, a_vp])
    plt.plot([0, b_vp], [0, 0])
    plt.plot([0, b_vp], [a_vp, 0])
    plt.plot([0, x_cm], [0, y_cm])

    plt.xlim(xmin=x_min, xmax=x_max)
    plt.ylim(ymin=y_min, ymax=y_max)

    plt.text(x_cm, y_cm, s=(str(round(cm, 1)) + r'$\degree$'), fontsize=16)
    plt.text(x_cm / 2, y_cm / 2, s=(str(round(v_cm, 1))))
    plt.text(x_cm, y_cm + n_radius, s=('R=' + str(round(n_radius, 1))))

    # First reactive circle
    circle = plt.Circle((x_cm, y_cm), radius=n_radius, fill=0, color='red')
    plt.gca().add_patch(circle)
    plt.text(x_cm + n_radius, y_cm, s=(str(round(high_limit, 1)) + r'$\degree$'), fontsize=16)
    plt.text(x_cm - n_radius, y_cm, s=(str(round(low_limit, 1)) + r'$\degree$'), fontsize=16)

    plt.text(x=(x_max - x_min) / 2, y=(y_max - y_min) / 2, s=(r'Ec = %s kJ/mol ' % round(ec / 1000, 1)), fontsize=12)
    plt.text(x=-100, y=a_vp / 2, rotation=90, s=(r'$v_{1{\degree}}$ = %s m/s ' % round(a_vp, 1)), fontsize=12)
    plt.text(x=b_vp / 2, y=-100, s=(r'$v_{2{\degree}}$ = %s m/s ' % round(b_vp, 1)), fontsize=12)
   
    #save figure as pdf 
    if save_pdf==True:
        plt.savefig('CM.pdf', format='pdf', bbox_inches='tight')

    # Let plot values to the convenience 
    plt.figure(figsize=(6, 6))
    plt.axis('off')
    plt.text(x=0.1, y=1.26, s=(f'{A}({a_mass} amu) + {B}({b_mass} amu) = {C}({c_mass} amu) + {D}({d_mass} amu)'),
             fontsize=26, weight='bold')
    plt.text(x=0.1, y=1.16, s=(r'$\Delta H\degree$ = - %s kJ/mol' % round(er / 1000, 2)), fontsize=16)

    plt.text(x=0.1, y=0.96, s=('CM angle = %s' % round(cm, 2) + r'${\degree}$'), fontsize=16, c='red')
    plt.text(x=0.1, y=0.86, s=('max = %s' % round(high_limit, 2) + r'${\degree}$'), fontsize=16)
    plt.text(x=0.1, y=0.76, s=('min = %s' % round(low_limit, 2) + r'${\degree}$'), fontsize=16)
    plt.text(x=0.1, y=0.66, s=(r'Ec = %s kJ/mol ' % round(ec / 1000, 1)), fontsize=16)
    plt.text(x=0.1, y=0.56,
             s=(r'$v_{1{\degree}}$ = %s m/s ' % round(a_vp, 1) + r' with $s_{1{\degree}}$ = %s' % round(a_s, 1)),
             fontsize=16)
    plt.text(x=0.1, y=0.46,
             s=(r'$v_{2{\degree}}$ = %s m/s ' % round(b_vp, 1) + r' with $s_{2{\degree}}$ = %s' % round(b_s, 1)),
             fontsize=16)
    plt.text(x=0.1, y=0.36, s=(r'$v_{cm}$ = %s m/s ' % round(v_cm, 1)), fontsize=16)
    plt.text(x=0.1, y=0.26, s=(r'TOF peak center = %s s ' % round(t_o_f, 1)), fontsize=16)
    plt.text(x=0.1, y=0.16, s=(r'TOF starts at = %s s ' % round(t_o_f_low, 1)), fontsize=16)
    plt.text(x=0.1, y=0.06, s=(r'TOF ends at = %s s ' % round(t_o_f_high, 1)), fontsize=16)
    
    print('peak center in channels = ' + str(t_o_f_ch))

def center_2ch(a_mass, b_mass, c_mass, d_mass, a_vp, a_s, b_vp, b_s,
               offset, er, A, B, C, D, A_carrier_gas,
               a1_mass, b1_mass, c1_mass, d1_mass, a1_vp, a1_s, b1_vp, b1_s, er_2, x_max=2200, x_min=-200, y_max=2200, y_min=-200, save_pdf=False):
    #This function solve CM angle, Ec and so on in case you need to overlay two CM cases together
    #a_mass, mass of a
     #b_mass,
     #c_mass,
     #d_mass,
     #a_vp, most prob vel of a
     #a_s, S of A beam
     #b_vp, most prob vel of B
     #b_s, S of B Beam
     #offset, offset to calculate TOF parameters
     #er, reaction energy
     #A, name of A
     #B, name of B
     #C, name of C
     #D, name of C
     #A_carrier_gas,
     #x_max = 2200, x on for CM fig
    #x_min = -200, x on for CM fig
    #y_max = 2200, y on for CM fig
    #y_min = -200 y on for CM fig
    # save_pdf=False = make True if you want to save a pdf with name CM_2Ch.pdf
    import math
    import matplotlib.pyplot as plt
    import numpy as np

    # For future LaTeX not to be ugly
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Arial'
    plt.rcParams['mathtext.it'] = 'Arial'
    plt.rcParams['mathtext.bf'] = 'Arial'
    plt.rc('font', family='arial', weight='normal', size=12)

    # offset of the detector
    cm_offset = 0.75

    # ec = energy of collision
    ec = 0.5 * a_mass * b_mass * 0.001 * (a_vp * a_vp + b_vp * b_vp) / (a_mass + b_mass)

    x = ((er + ec) * 2 * 1000 * (a_mass + b_mass) / c_mass / d_mass) ** 0.5
    n_radius = d_mass * x / (d_mass + c_mass)
    cm_speed = (a_mass ** 2 * a_vp ** 2 + b_mass ** 2 * b_vp ** 2) ** 0.5 / (a_mass + b_mass)

    print(n_radius)
    print(cm_speed)

    if n_radius / cm_speed >= 1.0:
        high_limit = math.atan(b_mass * b_vp / a_mass / a_vp) * 180 / math.pi + cm_offset + math.asin(1) * 180 / math.pi
        low_limit = math.atan(b_mass * b_vp / a_mass / a_vp) * 180 / math.pi + cm_offset - math.asin(1) * 180 / math.pi
    else:
        high_limit = math.atan(b_mass * b_vp / a_mass / a_vp) * 180 / math.pi + cm_offset + math.asin(
            n_radius / cm_speed) * 180 / math.pi
        low_limit = math.atan(b_mass * b_vp / a_mass / a_vp) * 180 / math.pi + cm_offset - math.asin(
            n_radius / cm_speed) * 180 / math.pi

    cm = math.atan(b_mass * b_vp / a_mass / a_vp) * 180 / math.pi + cm_offset
    # high_limit = math.atan(b_mass*b_vp/a_mass/a_vp)*180/math.pi+cm_offset+math.asin(n_radius/cm_speed)*180/math.pi
    # low_limit =  math.atan(b_mass*b_vp/a_mass/a_vp)*180/math.pi+cm_offset-math.asin(n_radius/cm_speed)*180/math.pi
    x_cm = a_vp / (math.tan(math.pi / 2 - math.atan(b_mass * b_vp / a_mass / a_vp)) + a_vp / b_vp)
    y_cm = x_cm * math.tan(math.pi / 2 - math.atan(b_mass * b_vp / a_mass / a_vp))
    # low_limit_x = a_vp/(math.tan(math.pi/2-math.atan(b_mass*b_vp/a_mass/a_vp)+math.asin(n_radius/cm_speed))+a_vp/b_vp)
    # low_limit_y = low_limit_x*math.tan(math.pi/2-math.atan(b_mass*b_vp/a_mass/a_vp)+math.asin(n_radius/cm_speed))
    # high_limit_x = a_vp
    # high_limit_y = high_limit_x*math.tan(math.pi/2-(math.atan(b_mass*b_vp/a_mass/a_vp)+math.asin(n_radius/cm_speed)))
    # high_limit_y

    # calculate v_cm of a product
    v_cm = math.sqrt(x_cm ** 2 + y_cm ** 2)

    # estimate time of flight of the product
    t_o_f = 329 / 1000 / v_cm * 10 ** 6 + offset
    t_o_f_ch = t_o_f / 10.24
    t_o_f_low = 329 / 1000 / (v_cm + n_radius) * 10 ** 6 + offset
    t_o_f_low_ch = t_o_f_low / 10.24
    t_o_f_high = 329 / 1000 / (v_cm - n_radius) * 10 ** 6 + offset
    t_o_f_high_ch = t_o_f_high / 10.24

    #calculate channel 2

    x_2 = ((er_2 + ec) * 2 * 1000 * (a1_mass + b1_mass) / c1_mass / d1_mass) ** 0.5
    n_radius_2 = d1_mass * x_2 / (d1_mass + c1_mass)

    # ec = energy of collision
    ec_2 = 0.5 * a1_mass * b1_mass * 0.001 * (a1_vp * a1_vp + b1_vp * b1_vp) / (a1_mass + b1_mass)

    x_2 = ((er_2 + ec) * 2 * 1000 * (a1_mass + b1_mass) / c1_mass / d1_mass) ** 0.5
    n_radius_2 = d1_mass * x_2 / (d1_mass + c1_mass)
    cm_speed_2 = (a1_mass ** 2 * a1_vp ** 2 + b1_mass ** 2 * b_vp ** 2) ** 0.5 / (a1_mass + b1_mass)

    print(n_radius_2)
    print(cm_speed_2)

    if n_radius_2 / cm_speed_2 >= 1.0:
        high_limit_2 = math.atan(b1_mass * b1_vp / a1_mass / a1_vp) * 180 / math.pi + cm_offset + math.asin(
            1) * 180 / math.pi
        low_limit_2 = math.atan(b1_mass * b1_vp / a1_mass / a1_vp) * 180 / math.pi + cm_offset - math.asin(
            1) * 180 / math.pi
    else:
        high_limit_2 = math.atan(b1_mass * b1_vp / a1_mass / a1_vp) * 180 / math.pi + cm_offset + math.asin(
            n_radius_2 / cm_speed_2) * 180 / math.pi
        low_limit_2 = math.atan(b1_mass * b1_vp / a1_mass / a1_vp) * 180 / math.pi + cm_offset - math.asin(
            n_radius_2 / cm_speed_2) * 180 / math.pi

    cm_2 = math.atan(b1_mass * b1_vp / a1_mass / a1_vp) * 180 / math.pi + cm_offset
    # high_limit = math.atan(b_mass*b_vp/a_mass/a_vp)*180/math.pi+cm_offset+math.asin(n_radius/cm_speed)*180/math.pi
    # low_limit =  math.atan(b_mass*b_vp/a_mass/a_vp)*180/math.pi+cm_offset-math.asin(n_radius/cm_speed)*180/math.pi
    x_cm_2 = a1_vp / (math.tan(math.pi / 2 - math.atan(b1_mass * b1_vp / a1_mass / a1_vp)) + a1_vp / b1_vp)
    y_cm_2 = x_cm_2 * math.tan(math.pi / 2 - math.atan(b1_mass * b1_vp / a1_mass / a1_vp))

    # now plot

    plt.figure(figsize=(10, 10))
    plt.plot([0, 0], [0, a_vp])
    plt.plot([0, b_vp], [0, 0])
    plt.plot([0, b_vp], [a_vp, 0])
    plt.plot([0, x_cm], [0, y_cm])

    plt.plot([0, 0], [0, a1_vp])
    plt.plot([0, b1_vp], [0, 0])
    plt.plot([0, b1_vp], [a1_vp, 0])
    plt.plot([0, x_cm_2], [0, y_cm_2])

    plt.xlim(xmin=x_min, xmax=x_max)
    plt.ylim(ymin=y_min, ymax=y_max)

    plt.text(x_cm, y_cm, s=(str(round(cm, 1)) + r'$\degree$'), fontsize=16)

    # First reactive circle
    circle = plt.Circle((x_cm, y_cm), radius=n_radius, fill=0, color='red')
    plt.gca().add_patch(circle)
    plt.text(x_cm + n_radius, y_cm, s=(str(round(high_limit, 1)) + r'$\degree$'), fontsize=12)
    plt.text(x_cm - n_radius, y_cm, s=(str(round(low_limit, 1)) + r'$\degree$'), fontsize=12)

    # first elastic circle
    circle2 = plt.Circle((x_cm_2, y_cm_2), radius=n_radius_2, fill=0, color='blue')
    plt.text(x_cm_2, y_cm_2, s=(str(round(cm_2, 1)) + r'$\degree$'), fontsize=16)
    plt.gca().add_patch(circle2)
    plt.text(x_cm_2 + n_radius_2, y_cm_2 + 50, s=(str(round(high_limit_2, 1)) + r'$\degree$'), fontsize=12)
    plt.text(x_cm_2 - n_radius_2, y_cm_2 + 50, s=(str(round(low_limit_2, 1)) + r'$\degree$'), fontsize=12)

    #plt.text(100, 100, s=('V_1 = 2200, V_2 = '), fontsize=12)

    #save figure as pdf 
    if save_pdf==True:
        plt.savefig('CM_2Ch.pdf', format='pdf', bbox_inches='tight')
                   
    plt.show()


def cm_dev(a_mass, b_mass, c_mass, d_mass, a_vp, d_a_vp, a_s, d_a_s, b_vp, d_b_vp, b_s, d_b_s):
    import statistics as st
    import math
    # this function in use when you need to calculate all possible deviations of Ec, V, CM angles 
    
    #A + B = C + D
    #a_mass, b_mass, c_mass, d_mass - molecular masses of A,B,C,D
    # a_vp, b_vp = Vp of A and B in m/s,
    # d_a_vp, d_b_vp =  std dev of A and B Vp
    # a_s, b_s = S parameter of for the A and B beams
    # d_a_s, d_b_s = Std dev for S parameters of for the A and B beams

    #calculate different Collisional Energies
    ec = 0.5 * a_mass * b_mass * 0.001 * (a_vp * a_vp + b_vp * b_vp) / (a_mass + b_mass)
    ec_high = 0.5 * a_mass * b_mass * 0.001 * ((a_vp + d_a_vp) ** 2 + (b_vp + d_b_vp) ** 2) / (a_mass + b_mass)
    ec_low = 0.5 * a_mass * b_mass * 0.001 * ((a_vp - d_a_vp) ** 2 + (b_vp - d_b_vp) ** 2) / (a_mass + b_mass)

    #create a list of diff Ec and find mean and st dev
    list_of_e_coll = [ec_low, ec, ec_high]
    ec_mean = st.mean(list_of_e_coll)
    std_ec = st.stdev(list_of_e_coll)

    #calculate all possible CM angles
    cm = math.atan(b_mass * b_vp / a_mass / a_vp) * 180 / math.pi
    cm_1 = math.atan(b_mass * (b_vp + d_b_vp) / a_mass / (a_vp + d_a_vp)) * 180 / math.pi
    cm_2 = math.atan(b_mass * (b_vp + d_b_vp) / a_mass / (a_vp - d_a_vp)) * 180 / math.pi
    cm_3 = math.atan(b_mass * (b_vp - d_b_vp) / a_mass / (a_vp + d_a_vp)) * 180 / math.pi
    cm_4 = math.atan(b_mass * (b_vp - d_b_vp) / a_mass / (a_vp - d_a_vp)) * 180 / math.pi

    #create a list of CM angles
    cm_list = [cm, cm_1, cm_2, cm_3, cm_4]
    #find mean and st dev based on lowest and highest possible CM
    cm_mean = st.mean([max(cm_list), min(cm_list)])
    std_cm = st.stdev([max(cm_list), min(cm_list)])

    #create a fun to round all items in the list and return the new list
    def round_list(list, significant_figures, ratio=1):
        num = 0
        list_2 = []
        for i in list:
            list_2.append(round(i/ratio, significant_figures))
            num += 1
        return list_2

    print(f'CM = {round(cm_mean, 2)} ± {round(std_cm, 2)}')
    print(f'Ec = {round(ec_mean / 1000, 2)} ± {round(std_ec / 1000, 2)} kJ/mol')
    print(f'CM list {round_list(cm_list, 2)}')
    print(f'Ec list {round_list(list_of_e_coll, 2, ratio=1000)}')
    print(f'Vp(A) = {a_vp} ± {d_a_vp} m/s')
    print(f'S(A) = {a_s} ± {d_a_s}')
    print(f'Vp(B) = {b_vp} ± {d_b_vp} m/s')
    print(f'S(B) = {b_s} ± {d_b_s}')


# TOF 
# Here starts the section for faster work and process with TOFs (time-of-flight) recorded during Xbeams experiments
# files must be in the same folder as your Python/Jupyter project
# Main description of parameters:
# laser_off=True or False - to show or now laser-off line
# x_max = 200, x channels to show
# 'all_raw' -  all raw files in a folder
# 'all_man' - all manipulated (laser-on minus laser-off part of TOF) files a folder
# 'AC66611.TOF'- plot 1 files
# ['AC66611.TOF', 'AC66711.TOF'] - list of files
#Examples of use:
# read_tof("AC762810.TOF", laser_off=True) #only one  file
# read_tof(["AC762810.TOF","AC762910.TOF"], laser_off=True) #all from list
# read_tof('all_raw', laser_off=True) #all from list
###############

import pandas as pd
from IPython.display import display, HTML
import numpy as np
import matplotlib.pyplot as plt


def get_the_shift(file_name):
    with open(file_name, "r") as file:
        line_number = 12
        # get the line
        shift = file.readlines()[line_number - 1]
        # remove 'Shift= ' from line
        shift = shift.removeprefix('Shift= ')
        # make int from str
        shift = int(shift)
    return shift


# define fun that will join all files
def join_tofs(raw_or_man):
    import pandas as pd
    import os
    import glob
    import openpyxl

    if raw_or_man == 'man':
        # find the current working directory
        path = os.getcwd()

        # make a list of all manipulated TOFs
        tof_files = glob.glob(os.path.join(path, "*11.TOF"))

        # create a dummy df that have indexes like first TOF file
        one_tof = pd.read_table(tof_files[0], header=0, skiprows=11, names=[tof_files[0]])
        all_man_tofs = pd.DataFrame(index=one_tof.index)

        # loop for all TOFs in the list
        for tof_1 in tof_files:
            short_name = tof_1.removeprefix(path).removeprefix('\\').removesuffix(
                '.TOF')  # short name for the df table header
            one_tof = pd.read_table(tof_1, header=0, skiprows=11, names=[short_name])  # read a file
            all_man_tofs = all_man_tofs.join(one_tof)  # join the readed file to a data frame

        return all_man_tofs
        # all_man_tofs.to_csv(os.getcwd()+'\\all_man_data.csv', index=False)
        # all_man_tofs.to_excel(os.getcwd()+'\\all_man_data.xlsx', index=False)

    # Raw Data - all the same as above only diff is the pattern of the file name - 10.TOF

    elif raw_or_man == 'raw':

        path = os.getcwd()
        tof_files = glob.glob(os.path.join(path, "*10.TOF"))
        one_tof = pd.read_table(tof_files[0], header=0, skiprows=11, names=[tof_files[0]])
        all_tofs = pd.DataFrame(index=one_tof.index)

        for tof in tof_files:
            short_name = tof.removeprefix(path).removeprefix('\\').removesuffix('.TOF')
            one_tof = pd.read_table(tof, header=0, skiprows=11, names=[short_name])
            all_tofs = all_tofs.join(one_tof)

        return all_tofs, tof_files


def read_tof(file_name, laser_off=True, skip=True, x_max=200):
    # laser_off: if true will make an additional column with laser off_signal,. Shift will be taken from a file
    # skip = if false will read a file without skipping lines - which is make sense onle if you deleted firs 11 lines y your hand

    if file_name == 'all_raw':
        tof, tof_files = join_tofs('raw')
        tof['Average_TOF'] = tof.mean(axis=1)

        shift = get_the_shift(tof_files[0])

        tof['Laser_off'] = pd.DataFrame(tof['Average_TOF'][shift - 1:].reset_index().drop(['index'], axis=1))

        plt.figure(figsize=(9, 6))
        plt.plot(tof['Average_TOF'], label='Average_TOF')
        plt.xlim(0, x_max)
        if laser_off == True:
            plt.plot(tof['Laser_off'], label='Laser_off')
        plt.legend()

        return tof

    elif file_name == 'all_man':
        tof = join_tofs('man')
        tof['Average_TOF'] = tof.mean(axis=1)

        plt.figure(figsize=(9, 6))
        plt.plot(tof['Average_TOF'], label='Average_TOF')
        plt.xlim(0, x_max)
        plt.legend()

        return tof


    elif type(file_name) == str:

        if skip == True:
            skiprows_num = 11
        else:
            skiprows_num = 0

        short_name = file_name.removesuffix('.TOF')

        tof = pd.read_table(file_name, header=0, skiprows=skiprows_num, names=[short_name])

        plt.figure(figsize=(9, 6))
        # plot laser_on
        plt.plot(tof.iloc[0:, 0], label=tof.columns[0])
        plt.xlim(0, x_max)

        if laser_off == True:
            # get the shift from a file
            with open(file_name, "r") as file:
                line_number = 12
                # get the line
                shift = file.readlines()[line_number - 1]
                # remove 'Shift= ' from line
                shift = shift.removeprefix('Shift= ')
                # make int from str
                shift = int(shift)
            # write a laser off col in a different
            tof['Laser_off'] = pd.DataFrame(tof.iloc[shift - 1:, 0]).reset_index().drop(['index'], axis=1)
            # plot laser_off signal
            plt.plot(tof.iloc[0:, 1], label=tof.columns[1])

        plt.legend()

        return tof

    elif type(file_name) == list:
        tof = pd.read_table(file_name[0], header=0, skiprows=11, names=[file_name[0]])
        for single_tof in file_name:
            if skip == True:
                skiprows_num = 11
            else:
                skiprows_num = 0

            short_name = single_tof.removesuffix('.TOF')

            single_tof_df = pd.read_table(single_tof, header=0, skiprows=skiprows_num, names=[short_name])

            tof = tof.join(single_tof_df)

        tof = tof.drop(file_name[0], axis=1)

        tof['Average_TOF'] = tof.mean(axis=1)

        plt.figure(figsize=(9, 6))
        plt.plot(tof['Average_TOF'], label='Average_TOF')
        plt.xlim(0, x_max)

        if laser_off == True:
            # get the shift from a file
            with open(file_name[0], "r") as file:
                line_number = 12
                # get the line
                shift = file.readlines()[line_number - 1]
                # remove 'Shift= ' from line
                shift = shift.removeprefix('Shift= ')
                # make int from str
                shift = int(shift)
            # write a laser off col in a different
            tof['Laser_off'] = pd.DataFrame(tof['Average_TOF'][shift - 1:].reset_index().drop(['index'], axis=1))
            # plot laser_off signal
            plt.plot(tof['Laser_off'], label='Laser_off')
        plt.legend()

        return tof


def get_noise_level_of_average(data_from_read_tof, noise_starts=1629, noise_ends=2048):
    #import is output (dataframe) from read_tof function that have 'Average_TOF' column
    #it will return a value of noise in chossed region 
    #find col_index of Avergae data
    col_index=data_from_read_tof.columns.get_loc('Average_TOF')
    #cut noise level region
    #you can change 
    noise_level = data_from_read_tof.iloc[noise_starts:noise_ends, col_index].values
    #average it region to get a noise level
    noise_level = np.sum(noise_level) / len(noise_level)
    return noise_level


def get_noise_level_of_column(column_link, noise_starts=1629, noise_ends=2048):
    #import is any column
    #it will return a value of noise in chossed region 
    #extract region of noise
    noise_level = column_link[noise_starts:noise_ends].values
    #average it region to get a noise level
    noise_level = np.sum(noise_level) / len(noise_level)
    return noise_level

#function add new column with Time data to dataframe
#since our TOFs recorded with step (channels) you need to know wthat was a step size of ezch channel. in 99% cases for reactive scattering it was 10.24 us
def add_tofs_column(df,step=10.24, col_name='TOF'):
    import numpy as np 
    df[col_name] = np.arange(0, len(df)*step, step)


#increase dwell time in the file by sum of every two recorded dwell periods in a single tof
#sum every two rows = make 2.56 dwell as 5.12
def dwell_x_2(file):
    tof = read_tof(file, x_max=300)
    #drop average column
    tof = tof.drop('Laser_off', axis=1)
    N = 2
    tof = tof.groupby(tof.index // N).sum()
    return tof

#this function returns convertion of v(t) to v(p) for xbeams. S is a speed ratio. See [1] and [2] at the beginning of the file for a better explanation of v(t), v(p), S 
def vt_to_vp(vt, s):
    from math import sqrt
    vp=vt*(s+ sqrt(s*s+4))/(s+ sqrt(s*s+8))
    return vp

#tthis function returns conversion of v(p) to v(t) for xbeams. S is a ratio. See [1] and [2] in the beginning of the file better explanation of v(t), v(p), S 
def vp_to_vt(vp, s):
    from math import sqrt
    vt=vp/((s+sqrt(s*s+4))/(s+sqrt(s*s+8)))
    return vt
