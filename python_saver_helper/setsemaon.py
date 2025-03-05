#!python

import scipy.io as sio
import os
from os.path import dirname, join as pjoin
import time
from datetime import datetime

def setsemaon(subject, session, step_no):

  print('setting semaphore on')
  now = datetime.now()
  print("now =", now)
  dt_string = subject + now.strftime("%Y%m%d_%H%M%S")

  music_data = os.getenv('MUSIC_DATA') 
  data_dir = pjoin(music_data, 'tmp')
  mat_fname_in = pjoin(data_dir, 'recordDataFile_prototype.mat')
  mat_fname_out = pjoin(data_dir, 'recordDataFile.mat')

  mat_contents = sio.loadmat(mat_fname_in)
  mat_contents['recordData']['stepType'][0][0] = 'd'
  mat_contents['recordData']['recordState'][0][0] = 1
  mat_contents['recordData']['subjectName'][0][0] = dt_string
  mat_contents['recordData']['sessionName'][0][0] = session
  mat_contents['recordData']['stepNo'][0][0] = step_no
  mat_contents['recordData']['stepSegNo'][0][0] = 1 
  # If we turn this off, sequence is reprogrammed to SD only mode
  # and changes for sdlarge are not compatible
  mat_contents['recordData']['UIModePriority'][0][0] = 1

  #       str(mat_contents['recordData']['UIModePriority'][0][0]))
  sio.savemat(mat_fname_out, mat_contents)

  filename_pre = session + '/' + dt_string
#  print("step type: " + str(mat_contents['recordData']['stepType'][0][0]) +


  return filename_pre

def main():
  setsemaon('', '')

if __name__ == "__main__":
    main()
