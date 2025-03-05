function writeh5vlstringfn(fid, h5Path, strings)

% proc.args, vl string
VLstr_type = H5T.copy('H5T_C_S1');
H5T.set_size(VLstr_type,'H5T_VARIABLE');
% Create a dataspace for cellstr
H5S_UNLIMITED = H5ML.get_constant_value('H5S_UNLIMITED');
dspace = H5S.create_simple(1,numel(strings),H5S_UNLIMITED);
% Create a dataset plist for chunking
plist = H5P.create('H5P_DATASET_CREATE');
H5P.set_chunk(plist,1); % 1 string per chunk
% Create dataset
f=1;
dset = H5D.create(fid, h5Path, VLstr_type, dspace, plist);
% Write data
H5D.write(dset, VLstr_type,'H5S_ALL','H5S_ALL','H5P_DEFAULT', ...
          strings);
% Close file & resources
H5D.close(dset);
H5P.close(plist);
H5T.close(VLstr_type);
H5S.close(dspace);