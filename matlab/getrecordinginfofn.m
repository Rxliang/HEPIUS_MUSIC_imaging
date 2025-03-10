function [successFlag] = getrecordinginfofn(semaphoreOpenedState)
    global semaphoreKey;
    global recordDataFileName;
    global recordData; 
              
    if ~isempty(semaphoreOpenedState) && (1 == semaphoreOpenedState)
        %Read data from shared file
        semaphore('wait',semaphoreKey);
        successFlag = 1;
        if exist(recordDataFileName, 'file')
          try
            tmp = load(recordDataFileName);
          catch
            successFlag = 0;
            disp([mfile ': sync error loading recordDataFileName. Probably being written right now. Implement python semaphore control. Wait for next time.'])
            return
          end
          
          recordData = tmp.recordData;
        else
          disp(['No recording file exists. Will not record IQ ' ...
                'data.']);
          recordData.sessionName = ''; 
          recordData.subjectName = '';           
          recordData.recordState = -1;	  
          recordData.stepNo = '';          
          recordData.stepType = '';                    
          recordData.trialNo = '';
          recordData.UIPriorityMode = 1;          
          save(recordDataFileName, 'recordData');
        end
        semaphore('post',semaphoreKey);
    else
        %disp(['Semaphore not opened. Not recording']);
        successFlag =  0;
    end
end

