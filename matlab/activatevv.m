if isdeployed
  return
end

pth = pwd;

[dum, osType] = unix('echo $OSTYPE');

cd(veraPath)
pwd
activate
vtHome=getenv('MUSIC_HOME');
matlabCodePath = fullfile(vtHome, 'matlab');
addpath(matlabCodePath);
cd(matlabCodePath);

