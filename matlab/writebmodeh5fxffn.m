function h5File = writebmodeh5fxffn(spec, proj, varName, outFilePrefix, ...
                                    type, instanceNumber, firstRun)

persistent cnt 

if nargin < 6
  firstRun = 1;
end

% frame by frame version
if nargin < 5
  type = 'single';
end

h5File = [outFilePrefix '.h5'];

if firstRun
  % write livetimeFrac
  fcpl = H5P.create('H5P_FILE_CREATE');
  fapl = H5P.create('H5P_FILE_ACCESS');
  fid = H5F.create(h5File, 'H5F_ACC_TRUNC', fcpl, fapl);

  plistg = 'H5P_DEFAULT';

  h5Path = '/header';
  gid = H5G.create(fid, h5Path, plistg,plistg,plistg);

  hdf5write(h5File, '/header/PDelta', spec.PData.PDelta, 'WriteMode', ...
          'append');

  hdf5write(h5File, '/header/Size', spec.PData.Size,  ...
            'WriteMode', 'append');
  
  hdf5write(h5File, '/header/Origin', spec.PData.Origin, ...
            'WriteMode', 'append');            
  
  hdf5write(h5File, '/header/lambda_mm', spec.lambda_mm, ...
            'WriteMode', 'append');            

  dateStr = datestr(now, 'yyyymmdd_HHMM');
  writeh5vlstringfn(fid, '/header/creationTime', cellstr(dateStr));
  H5G.close(gid);
  h5create(h5File, ['/data/' varName], [size(proj,1) size(proj,2) ...
                      Inf], 'ChunkSize', [size(proj,1) size(proj,2) ...
                      1]);
  h5create(h5File, ['/data/datenumFrame'], [Inf], 'ChunkSize', [1]);
  cnt(instanceNumber) = 1;
else
  h5write(h5File, ['/data/' varName], feval(type, proj), ...
            [1 1 cnt(instanceNumber)], ...
            [size(proj,1) size(proj,2) 1]);
  h5write(h5File, ['/data/datenumFrame'], spec.datenumFrame, ...
            cnt(instanceNumber), 1);      
  cnt(instanceNumber) = cnt(instanceNumber)+1;
  [instanceNumber cnt]
end

end

