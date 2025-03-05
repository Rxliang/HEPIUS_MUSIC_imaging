if isempty(veraPath) 
  error('Environment variable VANTAGE_HOME not defined on system.');
else
  if ~exist(veraPath, 'dir')
    error(['Directory ' veraPath ' does not exist on this system.']);
  end
end

if isempty(vtHome) 
  error('Environment variable MUSIC_HOME not defined on system.');
else
  if ~exist(veraPath, 'dir')
    error(['Directory ' vtHome ' does not exist on this system.']);
  end
end

if isempty(vtData) 
  error('Environment variable MUSIC_DATA not defined on system.');
  else
  if ~exist(vtData, 'dir')
    error(['Directory ' vtData ' does not exist on this system.']);
  end
end

