import numpy as np
import scipy
import scipy.io
import scipy.interpolate
import sys
import pylab as plt
from numpy.linalg import inv
import os
import matplotlib.pyplot as plt

print('An array of the intermediate frames you want to generate, for example\n'\
      'if you start with 40 frames and intend to generate 660 frames\n'\
      'generating a series of intermediate frames, you can include an array\n'\
      'as follows\n')

print('[40,80,165,330,660]\n')

print('which generates 40 frames in the first iteration, 80 in the second\n'\
      '165 in the third, 330 in the fourth and 660 in the last one\n')

arr_str=input('Include string with the frames array\n')

arr=np.array(arr_str.split('[')[1].split(']')[0].split(',')).astype(int)

k_low=arr[:-1]
k_up=arr[1:]

print('Provide directory to the file containing the first frames\n'\
    'indicated in the first position of the frames array\n'\
    'if your first position is 40, the folder should contain at least 40\n'\
    'frames\n')

dir_path=input('Indicate directory path\n')

print('All frames images should contain the same preffix and extension\n'\
        'with only the number of the frame varying\n')

pref=input('Provide the preffix of the images frames\n')
preff_gen=dir_path+pref

print('##########Execution started############')
for n in np.arange(0,len(k_low),1):
    if n==0:
        lista=np.arange(1,k_low[n]+1,1)
        data={}
        image_fin={}

        for i in lista:
            data[i]=scipy.io.loadmat(preff_gen+str(i)+'.mat')
            for k in data[i].keys():
                if type(data[i][k])==np.ndarray:
                    key_arr=k
            image_fin[i]=data[i][key_arr]

    else:
        image_fin=im_fin
        lista=np.arange(1,k_low[n]+1,1)

    print('#########Interpolating frames data#########')
    vec_im=np.linspace(1,k_low[n],k_up[n])
    shape=np.shape(image_fin[1])
    im_fin={}
    for k in np.arange(1,len(vec_im)+1,1):
        im_fin[k]=np.zeros(shape)
    for i in np.arange(0,shape[0],1):
        for j in np.arange(0,shape[1],1):
            inter=np.array([])
            for k in lista:
                inter=np.append(inter,image_fin[k][i,j])
                inter[np.isnan(inter)]=0
            int_im=scipy.interpolate.interp1d(lista,inter)
            for l in np.arange(0,len(vec_im),1):
                im_fin[l+1][i,j]=int_im(vec_im[l])

    print('###########Saving interpolated frames########')
    dir_path=dir_path.split('/')[0]+'/frames_'+str(k_up[n])+'/'
    preff_gen=dir_path+pref
    os.system('mkdir '+dir_path)


    for i in np.arange(1,len(im_fin)+1,1):
        mdic={"mat":im_fin[i],"label":'mat'}
        scipy.io.savemat(preff_gen+str(i)+'.mat',mdic)

print('############Execution ended################')
