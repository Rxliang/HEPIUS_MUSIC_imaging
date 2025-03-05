mfile = mfilename;
vtData = getenv('MUSIC_DATA');

folder_list = {'0815Survival', '815Survival', '816dayone_survival', ...
               '816dayonesurvival1PM', '816dayonesurvival5PM', ...
               '816dayonesurvival9PM', '817daytwosuvival9AM', ...
               '817daytwosurvival9PM', '818daythreesurvival9AM', ...
               '818daythreesurvival9PM', '819dayfoursurvival12PM', ...
               '1107closeloopbaseline', '1107closeloopfus1', '1107fus3', ...
               '1107fusrep', '1107fussteer', '1219baseline', ...
               '1219_PostFusSession1_ClosedLoop', '1219_PostFusSession2_ClosedLoop', ...
               '1219_PostFUSSession3_ClosedLoop', '1219_EpiduralVessel_9PM'};

%folder_list = {'test'};

set_names = [];
sn=0;

for q = 1:length(folder_list)
    d = dir(fullfile(vtData, folder_list{q}, '*sg1.re'));
    for r = 1:length(d)
        set_name = split(d(r).name, 'sdl_s');
        sn=sn+1;
        spl =  split(set_name{2}, '_sg');
        set_name_suff = spl{1};
        
        set_names{sn}.name = [set_name{1} 'sdl_s' set_name_suff];                             
        dr = dir(fullfile(vtData, folder_list{q}, [set_names{sn}.name '_*sg*.re']));                
        set_names{sn}.num_seg = length(dr);
        set_names{sn}.folder = folder_list{q};
               
        for s = 1:set_names{sn}.num_seg
          fn = fullfile(vtData, folder_list{q}, [set_names{sn}.name '_sg' num2str(s) '.im']);
          if ~exist(fn, 'file')
            error(['Missing file: ' fn]);            
          else
              dir(fn);
          end
          
        end % segs
        
    end % r unique sets in folder
end % q folders

music_sets = [];
music_sets.set_names = set_names;
music_sets.folder_list = folder_list;
music_sets.root_folder = vtData;

dateStr = datestr(now, 'yyyymmdd');

mat_file = fullfile(vtData, [mfile '_' dateStr '.mat']);
save(mat_file, 'music_sets');
lslrt(mat_file);


