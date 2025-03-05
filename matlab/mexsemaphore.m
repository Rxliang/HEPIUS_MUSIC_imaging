pwdSave = pwd;
mPath = [getenv('BNS_HOME') '\matlab'];
cd(mPath);

mex semaphore.c
cd(pwdSave);



              
  
