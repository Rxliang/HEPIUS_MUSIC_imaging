
import scipy.io as sio
import os
from os.path import dirname, join as pjoin
import time

music_data = os.getenv('MUSIC_DATA') 

data_dir = pjoin(music_data, 'tmp')

mat_fname = pjoin(data_dir, 'recordDataFile.mat')

while 1:
  mat_contents = sio.loadmat(mat_fname)
  time.sleep(1)
  print("subjectName: " + str(mat_contents['recordData']['subjectName'][0][0]) +
        " step type: " + str(mat_contents['recordData']['stepType'][0][0]) +
        "UIM" + str(mat_contents['recordData']['UIModePriority'][0][0]) +
        "recordState" + str(mat_contents['recordData']['recordState'][0][0]))  
        
  
  
