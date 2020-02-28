import vtk
import os
import numpy as np
import matplotlib.pyplot as plt 
import cv2
import time
import shutil
import imageio


#####################################################
def load_data(folder_list):
    print("Loading vtp files ...")
    all_data = []
    for item in folder_list:
        total_data = []
        for filename in sorted(os.listdir(item)):
            if "particles_" in filename and ".vtp" in filename:
                reader = vtk.vtkXMLPolyDataReader()
                reader.SetFileName(os.path.join(item, filename))
                reader.Update()
                pdata = reader.GetOutput()
                total_data.append(pdata)
        all_data.append(total_data)

    print("Done.")

    return all_data


#####################################################
def get_fig(num,total_data, save_folder=None):
    v_data = total_data[num].GetPointData().GetArray('Velocity')
    direction = True
    v_x = []
    v_y = []
    x = []
    y = []
    num_particle = total_data[num].GetNumberOfPoints()
    for i in range(num_particle):
        v_x.append(v_data.GetTuple(i)[0])
        v_y.append(v_data.GetTuple(i)[1])
        x.append(total_data[num].GetPoint(i)[0])
        y.append(total_data[num].GetPoint(i)[1])
    count_negtive = 0
    count_positive = 0
    count_zeros = 0
    for i in range(num_particle):
        if v_x[i] > 0:
            count_positive += 1
        elif v_x[i] == 0:
            count_zeros += 1
        else:
            count_negtive += 1
    if (count_positive)/num_particle >= (count_negtive)/num_particle:
        direction = True
    else:
        direction = False
    v_x = np.array(v_x)
    v_y = np.array(v_y)
    x = np.array(x)
    y = np.array(y)
    bottom = int(num_particle*0.005)
    mid = int(bottom/2) 
    temp_x = np.array([])
    temp_y = np.array([])
    # for the first 10 iterations zeros
    
    v_total = np.sqrt(v_x**2 + v_y**2)
    
    if count_zeros/num_particle >= 0.6 and direction == True:
        y_peak_position=(y[v_total == v_total.max()].mean())
        x_peak_position=(x[v_total == v_total.max()].mean())
        max_v = 0
    else:
        if direction == True:
            for i in range(num_particle):
                if v_x[i] < 0:
                    v_x[i] = 0
                    v_y[i] = 0
        else:
            for i in range(num_particle):
                if v_x[i] > 0:
                    v_x[i] = 0
                    v_y[i] = 0
        v_total = np.sqrt(v_x**2 + v_y**2)
        
        temp_x = np.argsort(-v_total)[:bottom]
        temp_x = temp_x[mid-10:mid+10]
        temp_y = np.argsort(-v_total)[:bottom]
        temp_y = temp_y[mid-10:mid+10]
        x_peak_position = x[temp_x].mean()
        y_peak_position = y[temp_y].mean()
        max_v = v_total[temp_x]
    
    if save_folder is not None:
        fig = plt.figure() 
        ax1 = fig.add_subplot(111) 
        ax1.set_title('Scatter Plot') 
        plt.xlabel('X') 
        plt.ylabel('Y') 
        ax1.scatter(x,y) 
        ax1.scatter(x_peak_position,y_peak_position,s=128,c='r',marker='o') 
        plt.legend('x1') 
        fig1 = plt.gcf()
        fig1.savefig(os.path.join(save_folder, 'fig_{:05}.png'.format(num)), dpi=100)
        plt.close(fig1)

    return max_v,x_peak_position,y_peak_position


#####################################################
def get_x_y_v(total_data, save_folder=None):
    max_v = []
    max_x = []
    max_y = []
    total_information = []
    for i in range(len(total_data)):
        v,x,y = get_fig(i,total_data, save_folder)
        max_v.append(v)
        max_x.append(x)
        max_y.append(y)
    
    total_information.append(max_x)
    total_information.append(max_y)
    total_information.append(max_v)
    
    return total_information


#####################################################
def make_animation(figures_folder):
            
    with imageio.get_writer(os.path.join(figures_folder, 'animation.gif'), mode='I') as writer:
        for filename in sorted(os.listdir(figures_folder)):
            if ".png" in filename:
                image = imageio.imread(os.path.join(figures_folder, filename))
                writer.append_data(image)
    

############################################################################
if __name__ == "__main__":

    this_dir = os.path.dirname(os.path.abspath(__file__))
    folders_names = ["rectangle_20s_dx=5e-2", "run_008", "rectangle_20s_1e-1", "run_016", "rectangle_20s_2e-1", "run_04","run_05", "run_08"]

    folder_list = [os.path.join(this_dir, "data", name) for name in folders_names]

    all_data = load_data(folder_list)

    labels = ['0.05', '0.08', '0.1', '0.16', '0.2', '0.4', '0.5', '0.8']

    figures_folder = os.path.join(this_dir, "data", "figures")
    if os.path.exists(figures_folder):
        shutil.rmtree(figures_folder)
    os.mkdir(figures_folder)

    save_folders = [os.path.join(figures_folder, label) for label in labels]
    for folder in save_folders:
        os.mkdir(folder)

    total = []
    for i in range(len(folders_names)):
        t = []
        for j in range(len(all_data[i])):
            t.append(j/10)
        information = get_x_y_v(all_data[i], save_folders[i])
        information = get_x_y_v(all_data[i])
        total.append(np.array(information[0]))
    
    t = np.array(t)
    for folder in save_folders:
        make_animation(folder)
    
    fig, ax = plt.subplots(4, 2, figsize=(10, 10))
    plt.tight_layout()
    for i, ax_i in enumerate(ax.ravel()):
        ax_i.plot(t, total[i], label=labels[i])
        ax_i.legend()
    
    plt.savefig(os.path.join(this_dir, "present", "crest_x_pos.pdf"))
    plt.close(fig)

    slope_measure_ranges = [(6, 24), (92, 110), (175, 195)]
    slopes = [[], [], []]
    for x_values in total:
        for i, (start, end) in enumerate(slope_measure_ranges):
            slopes[i].append(np.polyfit(t[start:end], x_values[start:end], 1)[0])
    
    fig = plt.figure()
    for i, s in enumerate(slopes):
        plt.plot([float(l) for l in labels], s, label="{:d}. pass".format(i+1));
    plt.xlabel('mesh size')
    plt.ylabel('wave speed [m/s]')
    plt.xscale('log')
    plt.legend()
    plt.savefig(os.path.join(this_dir, "present", "crest_speed.pdf"))
