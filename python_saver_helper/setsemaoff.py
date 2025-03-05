import scipy.io as sio
import os
from os.path import dirname, join as pjoin
import time
from datetime import datetime

def setsemaoff():
  print('setting semaphore off')
  dt_string = ''
  music_data = os.getenv('MUSIC_DATA') 
  data_dir = pjoin(music_data, 'tmp')
  #mat_fname_in = pjoin(data_dir, 'recordDataFile_prototype.mat')
  mat_fname_in = pjoin(data_dir, 'recordDataFile.mat')
  mat_fname_out = pjoin(data_dir, 'recordDataFile.mat')
  mat_contents = sio.loadmat(mat_fname_in)
  print("session: " + str(mat_contents['recordData']['sessionName'][0][0]))
  print("subject: " + str(mat_contents['recordData']['subjectName'][0][0]))
  mat_contents['recordData']['stepType'][0][0] = ' '
  mat_contents['recordData']['recordState'][0][0] = 0
  # need to change the step number back to zero to trigger a final save
  # of the buffer
  mat_contents['recordData']['stepNo'][0][0] = 0
  mat_contents['recordData']['stepSegNo'][0][0] = 0  
#  mat_contents['recordData']['subjectName'][0][0] = dt_string
  
  sio.savemat(mat_fname_out, mat_contents)

#  print("step type: " + str(mat_contents['recordData']['stepType'][0][0]) +
 #       str(mat_contents['recordData']['UIModePriority'][0][0]))

 
def main():
  setsemaoff()

if __name__ == "__main__":
    main()
  
