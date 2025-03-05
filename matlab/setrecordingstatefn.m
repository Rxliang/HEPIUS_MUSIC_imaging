%%%%%%%%%%%%%%%%%%%%%% Helper function : setRecordingState  %%%%%%%%%%%%%%%%%%
% recordState : % -1 = Not initialized, 0 = Stopped 2 = Paused 1 = Started 3 = Resumed
function setrecordingstatefn(recordState)
global semaphoreKey;
global recordDataFileName;

recordData = evalin('base', 'recordData');

recordData.recordState = recordState; 
%write data in shared file
semaphore('wait',semaphoreKey);
save(recordDataFileName,'recordData');
semaphore('post',semaphoreKey);
end

