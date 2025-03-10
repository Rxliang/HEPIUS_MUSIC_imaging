if isdeployed
  return
end

vtHome=getenv('MUSIC_HOME');
matlabCodePath = fullfile(vtHome, 'matlab');
addpath(matlabCodePath);

hn = gethostnamefn;

if ~strcmp(hn, 'FULL-VVRILHIUJE')
  vantageDir = 'Vantage-4.9.2-2308102000';
else    
  vantageDir = 'Vantage-4.8.4-2211151000';
end

if strcmp(hn, 'FULL-HFMK7TQI8O')
    vantageDir = 'Vantage-4.9.2-2308102000';
    %  vantageDir = 'Vantage-4.5.3-2107301223';
end


if strcmp(hn, 'WINDOWS-RQO8E6B') 
    vantageDir = 'VantageNXT-1.0.4-p1';
    %    vantageDir = 'VantageNXT-1.1.0-p1';
end

if strcmp(hn, 'Verasonics-PC')  
    vantageDir = 'VantageNXT-1.1.0-p1';
end

if ~ispc
  [dum, osType] = unix('echo $OSTYPE');
  homeDir = getenv('HOME');
  srcRoot = homeDir;
  if strcmp(hn, 'gray')
    vantageDir = 'Vantage-4.9.2-2308102000';
    veraPath = fullfile(homeDir, 'verasonics', vantageDir);
  else
    vantageDir = 'VantageNXT-1.0.4-p1'
  end
else
  dum=0;
  osType = 'windows';
  homeDir = getenv('userprofile');
  veraPath = fullfile(homeDir, 'Documents', vantageDir);
  srcRoot = [homeDir '/Documents/'];              
end

if dum==0 & ~isdeployed
%  addpath([veraPath '/Tools/RFDataViewer']);
  addpath([veraPath '/Utilities/']);
  addpath([veraPath '/System/']);
  addpath([veraPath '/vstoolbox/']);  
%  addpath([veraPath '/Examples_Biomedical/Specialty_Applications/ColorDopplerImagingExternal/']);    
  addpath(veraPath);
end

if ~isdeployed
    cd(matlabCodePath);
end

activatevv


