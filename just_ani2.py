import shutil
from datetime import datetime
import os
import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
# import random
# from time import sleep


space = 64
num_procs = 4
total_run_time = 5.0

location = "/home/himanshu/Desktop/final_year_now/data/"


root_save_folder = datetime.now().strftime("%Y-%m-%d__%H:%M:%S")
os.mkdir(root_save_folder)
shutil.copyfile("./main.c", root_save_folder+"/main.c")
shutil.copyfile("./params.h", root_save_folder+"/params.h")


for subdir, dirs, files in os.walk(location):
    for folder in dirs:
        data_folder = os.path.join(location, folder)

        save_folder = root_save_folder+"/"+folder

        os.mkdir(save_folder)  # makes folder with names like run1

        mat_files = []
        mat_data = {}
        tse_files = []
        tse_data = {}

        for i in range(num_procs):
            mat_files.append(data_folder + "/"+"MAT" + str(i)+".dat")
            tse_files.append(data_folder + "/"+"TserE" +
                             "{:02}".format(i)+".dat")

        for i in range(num_procs):
            mat_data[i] = np.loadtxt(mat_files[i])
            tse_data[i] = np.loadtxt(tse_files[i])

        f, arr_plots = plt.subplots(num_procs, 1, figsize=(1, 1*num_procs))

        for i in range(num_procs):
            arr_plots[i].plot(tse_data[i])
            if i == 0:
                arr_plots[0].set_title("top")

        for ax in arr_plots.flat:
            ax.set(ylabel='Potential (nV)', xlabel='time (some_unit)')
            ax.label_outer()

        f.savefig('TSE.png')
        shutil.move("./TSE.png", save_folder+"/TSE.png")

        spacey = space
        cols_per_pro = (int)(space/num_procs)
        # print(cols_per_pro)

        number_of_frames = (int)(mat_data[0].shape[0]/cols_per_pro)
        # print(number_of_frames)
        frames = np.zeros((number_of_frames, space, spacey))

        for i in range(number_of_frames):

            frames[i] = np.vstack((mat_data[0][i*cols_per_pro:(i+1)*cols_per_pro, :],
                                   mat_data[1][i *
                                               cols_per_pro:(i+1)*cols_per_pro, :],
                                   mat_data[2][i *
                                               cols_per_pro:(i+1)*cols_per_pro, :],
                                   mat_data[3][i *
                                               cols_per_pro:(i+1)*cols_per_pro, :],
                                   ))

        fig = plt.figure()

        ims = []

        a = 2
        for i in range((int)(number_of_frames/a)):

            im = plt.imshow(frames[a*i], vmin=-55, vmax=30, animated=True)
            # print("Done ",i)
            if i == 0:
                plt.colorbar()
            ims.append([im])

        # ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
        #                                 repeat_delay=1000)

        ani = animation.ArtistAnimation(
            fig, ims, interval=100, blit=True, repeat=False)

        ani.save('video.mp4')
        shutil.move("./video.mp4", save_folder+"/video.mp4")
