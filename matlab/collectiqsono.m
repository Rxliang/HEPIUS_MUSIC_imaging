mfile = mfilename;
vtData = getenv('MUSIC_DATA');
% survival studies, 202408
%set_list_file = 'getmusiciqsets_20240903.mat';
%load(fullfile(vtData, set_list_file), 'music_sets');

%map_file = 'offlinesdlarge_matfilemap_getmusiciqsets_20241213.mat';
map_file = 'offlinesdlarge_matfilemap_getmusiciqsets_20241227.mat';
report_file = 'music_dataset_report_20241228.xlsx';

dateStr = datestr(now, 'yyyymmdd');

out_script = fullfile(vtData, [mfile '_' dateStr '.sh']);

load(fullfile(vtData, map_file), 'map');

out_path = ['gray:/down/iqsono'];
%rsync_bin = 'h:\Documents\cwrsync\bin\rsync';
cyg_bin = 'c:\cygwin64\bin';
rsync_bin = fullfile(cyg_bin, 'rsync');
rsync_bin = 'rsync';
exec_strs = [rsync_bin ' -azvP '];

for q = 1:length(map.sel)
  in_path = fileparts(map.matFile{map.sel(q)});
  [ret, cyg_path] = unix(fullfile(cyg_bin, ['cygpath ' in_path]));
  cyg_path = cyg_path(1:end-1); % ascii code 10 at end
  load(map.matFile{map.sel(q)}, 'topInd', 'xzRect')
   
  %'h:\Documents\cwrsync\bin\rsync'
  sono_file = [map.iq2sono{map.sel(q)} '.mat'];
  load(fullfile(in_path, sono_file), 'spec');
  num_img = length(spec.inFilesPre);
  im = [];
  for r = 1:1
    img_file = spec.inFilesPre{r};
    im_mat = readimfn(img_file, in_path, 1, 0);
    iset = load(im_mat);
    im(:,:,r) = iset.iSet.im;
  end

  figure(1)
  clf
  imagesc(iset.iSet.xAx_mm, iset.iSet.zAx_mm, im(:,:,1))
  hold on
  xzr = xzRect{map.sel(q)};
  hs = drawsquare(xzr(1,1), xzr(1,2), xzr(2,1), xzr(2,2));
  set(hs , 'color', 'r')
  ylim([2 12]);
  xlim([-10 10]);
  xlabel('x (mm)');
  ylabel('z (mm)');
  title(titstrfn(img_file))
  png_file = fullfile(in_path, [img_file '.png']);
  if 1 | ~exist(png_file, 'file')
    print('-dpng', '-f1', png_file);
    lslrt(png_file)
  end
  
  sono_path = [cyg_path '/' sono_file];
  img_path = [cyg_path '/' img_file '.png'];
  
  %  lslrt(sono_path)
  %  exec_str = [rsync_bin ' -azvP ' sono_path ' ' out_path '; ']
  exec_str = [sono_path ' ' img_path ' '];
  %  unix(exec_str)
  exec_strs = horzcat(exec_strs, exec_str);
  
end

report_path = fullfile(vtData, report_file);
[ret, cyg_path] = unix(fullfile(cyg_bin, ['cygpath ' report_path]));
cyg_path = cyg_path(1:end-1); % ascii code 10 at end
  
exec_strs = horzcat(exec_strs, [cyg_path ' ']);
exec_strs = horzcat(exec_strs, out_path);

%exec_str = [rsync_bin ' -azvP ' report_path ' ' out_path '; ']
%exec_strs = horzcat(exec_strs, exec_str);

fid = fopen(out_script, 'w');
fprintf(fid, '%s\n', exec_strs);
fclose(fid)
lslrt(out_script);

