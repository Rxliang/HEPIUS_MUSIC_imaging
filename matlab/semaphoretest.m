global semaphoreKey; semaphoreKey = 24537;
vtHome = getenv('VITTAMED_HOME');
global recordDataFileName; 
recordDataFileName = [vtHome '/matlab/tmp/recordDataFile.mat'];
global recordData;
global nonWindowsHost;
%persistent semaphoreOpenedState; %This has been changed from gloabal to persistent
                                 % each SD instance create / (mainly opens)
                                 % semaphore only once
%First time
semaphoreOpenedState = [];
if isempty(semaphoreOpenedState) || (0 == semaphoreOpenedState)
  [dum, osType] = unix('echo $OSTYPE');

  if dum==0 & (strcmp(osType(1:6), 'darwin') | strcmp(osType(1:5), ...
                                                      'linux'))
    nonWindowsHost=1;
  else
    nonWindowsHost=0;
    %Create OR OPEN semaphore - we use the mexfile associated with semaphore.c
    semaphore('create',semaphoreKey,1); %Initial semaphore value 1;
  end
  
  semaphoreOpenedState = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%

while 1
 if ~isempty(semaphoreOpenedState) && (1 == semaphoreOpenedState)
        %Read data from shared file
        semaphore('wait',semaphoreKey);
        successFlag = 1;
        if exist(recordDataFileName, 'file')
          tmp = load(recordDataFileName);
          recordData = tmp.recordData
        else
          disp(['No recording file exists. Will not record IQ ' ...
                'data.']);
          recordData.subjectName = ''; 
          recordData.recordState = -1;
          save(recordDataFileName, 'recordData');
        end
        semaphore('post',semaphoreKey);
    else
        disp(['Semaphore not opened. Not recording']);
        successFlag =  0;
 end
 pause(1)
end
