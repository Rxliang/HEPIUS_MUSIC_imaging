mfile = mfilename;
vtData = getenv('MUSIC_DATA');

folder_list = {'1219baseline/old'};
new_name = {'baseline'};

%folder_list = {'test'};

exec = 1;
set_names = [];
sn=0;

for q = 1:length(folder_list)
  dfix = dir(fullfile(vtData, folder_list{q}, '*_fix*.re'));
  fix_vec = [];
  for f = 1:length(dfix)
    tmp = split(dfix(f).name(1:end-3), '_fix');
    fix_vec(f) = str2num(tmp{2});
  end
  fix_vals = unique(fix_vec);

  for f = 1:length(fix_vals)
    d = dir(fullfile(vtData, folder_list{q}, ['*sg1_fix' num2str(fix_vals(f)) '.re']));
    for r = 1:length(d)
        set_name = split(d(r).name, 'sdl_s');
        sn=sn+1;
        spl =  split(set_name{2}, '_sg');
        set_name_suff = spl{2};        
        set_names{sn}.name = [set_name{1} 'sdl_s_sg' set_name_suff];
        dr = dir(fullfile(vtData, folder_list{q}, [set_name{1} 'sdl_s_sg*_fix' num2str(fix_vals(f)) '.re']));                
        set_names{sn}.num_seg = length(dr);
        set_names{sn}.folder = folder_list{q};
               
        for s = 1:set_names{sn}.num_seg
          fn_pre = fullfile(vtData, folder_list{q}, [set_name{1} 'sdl_s_sg' num2str(s) '_fix' num2str(fix_vals(f))]);
          fn_re = [fn_pre '.re'];
          fn_im = [fn_pre '.im'];
          fn_img = [fn_pre '.img'];          
          dst = fullfile(vtData, folder_list{q}, '..', [new_name{q} '_sdl_s' num2str(fix_vals(f)+1) '_sg'  num2str(s)]);

          disp(['original: ' fn_pre ' dst: ' dst])
          
          if ~exist(fn_re, 'file') | ~exist(fn_im, 'file')
            error(['Missing file: ' fn]);            
          else
            %dir(fn_re);
          end
          
          dst_re = [dst '.re'];
          dst_im = [dst '.im'];
          dst_img = [dst '.img'];
          
          if exist(dst_re, 'file') | exist(dst_im, 'file') | exist(dst_img, 'file')
            %error(['Destination file exists: ' dst]);
            disp(['Destination file exists: ' dst]);
          end
          if 1
            exec_str = ['copy ' fn_re ' ' dst_re ]
            if exec
              unix(exec_str)
            end
            exec_str = ['copy ' fn_im ' ' dst_im]
            if exec
              unix(exec_str)
            end
            exec_str = ['copy ' fn_img ' ' dst_img]
            if exec
              unix(exec_str)
            end            

            
          end
        end % segs        
    end % r unique sets in folder
  end % fix vals
  
end % q folders

  return
  
music_sets = [];
music_sets.set_names = set_names;
music_sets.folder_list = folder_list;
music_sets.root_folder = vtData;

dateStr = datestr(now, 'yyyymmdd');

mat_file = fullfile(vtData, [mfile '_' dateStr '.mat']);
save(mat_file, 'music_sets');
lslrt(mat_file);


