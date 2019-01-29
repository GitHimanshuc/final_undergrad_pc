import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os
from datetime import datetime
import shutil
# import random
# from time import sleep



space = 64
num_procs = 4
run_number = 1


location = "/home/himanshu/Desktop/final_year_now/data/"
run1 =location+"run1/TserE00.dat"
run2 =location+"run1/TserE01.dat"
run3 =location+"run1/TserE02.dat"
run4 =location+"run1/TserE03.dat"


mat_files = []
mat_data = {}


for i in range(num_procs):
    mat_files.append(location +"run" + str(run_number) + "/"+"MAT"+ str(i)+".dat")
    
for i in range(num_procs):
    mat_data[i] = np.loadtxt(mat_files[i])


data1 = np.loadtxt(run1)
data2 = np.loadtxt(run2)
data3 = np.loadtxt(run3)
data4 = np.loadtxt(run4)


f, arr_plots = plt.subplots(2,2)
arr_plots[0,0].plot(data1)
arr_plots[0,0].set_title("1 (top)")
arr_plots[0,1].plot(data2)
arr_plots[0,1].set_title("2")
arr_plots[1,0].plot(data3)
arr_plots[1,0].set_title("3")
arr_plots[1,1].plot(data4)
arr_plots[1,1].set_title("4 (bottom)")

for ax in arr_plots.flat:
    ax.set(ylabel = 'Potential (nV)', xlabel = 'time (some_unit)')
    ax.label_outer()



f.savefig('TSE.png')








spacey = space
cols_per_pro = (int)(space/num_procs)
print(cols_per_pro)

number_of_frames = (int)(mat_data[0].shape[0]/cols_per_pro)
print(number_of_frames)
frames = np.zeros((number_of_frames , space , spacey))

for i in range(number_of_frames):

    frames[i] = np.vstack((mat_data[0][i*cols_per_pro:(i+1)*cols_per_pro,:],
                           mat_data[1][i*cols_per_pro:(i+1)*cols_per_pro,:],
                           mat_data[2][i*cols_per_pro:(i+1)*cols_per_pro,:],
                           mat_data[3][i*cols_per_pro:(i+1)*cols_per_pro,:],
                          ))


fig = plt.figure()


ims = []

a = 2
for i in range((int)(number_of_frames/a)):

    im = plt.imshow(frames[a*i],vmin = -55,vmax = 30, animated=True)
    print("Done ",i)
    if i == 0:
      plt.colorbar()
    ims.append([im])



# ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
#                                 repeat_delay=1000)

ani = animation.ArtistAnimation(fig, ims, interval=100, blit=True,repeat=False)

ani.save('video.mp4')


a = datetime.now().strftime("%Y-%m-%d__%H:%M:%S")
os.mkdir(a)


shutil.copyfile("./main.c", a+"/main.c")
shutil.copyfile("./params.h",a+"/params.h")
shutil.move("./video.mp4",a+"/video.mp4")
shutil.move("./TSE.png",a+"/TSE.png")





