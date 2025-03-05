%global recordData

dateStr = datestr(now, 'yyyymmdd_HHMM');
recordData.subjectName = dateStr;
recordData.recordState = 1;
recordData.sessionName = 'semaphoretest';
recordData.stepNo = 1;
recordData.stepTrialNo = 1;
recordData.stepType = 'd';
recordData.UIModePriority = 1;

setrecordingstatefn(recordData.recordState)