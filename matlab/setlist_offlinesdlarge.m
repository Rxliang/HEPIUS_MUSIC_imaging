
if 0
  inPath = fullfile(vtData, 'sss');
  inFilePrePre = 'sss20240307_003351_sdl_s3_sg1';

  inPath = fullfile(vtData, 'pq');
  inFilePrePre = 'ps20240314_120621_sdl_s7_sg1';
  inFilePre = [inFilePrePre]; 
  r=10; c=10;
  
  inPath = fullfile(vtData, 'string');
  inFilePrePre = 'anon20240317_180125_sdl_s3_sg3';
  inFilePre = [inFilePrePre]; 
  r=1:10; c=1:5;

  inPath = fullfile(vtData, 'nxt');
  inFilePrePre = 'string5x1020240319_071929_sdl_s9_sg1';
  inFilePre = [inFilePrePre]; 
  r=1:10; c=1:5;
  r=5; c=4:5;
  
  inPath = fullfile(vtData, 'nxt_new_music');
  inFilePrePre = 'radial20x1020240321_141005_sdl_s11_sg1';
  inFilePre = [inFilePrePre]; 
  %r=1:10; c=1:5;
  %r=10:11; c=5;
  
  inPath = fullfile(vtData, 'string');
  inFilePrePre = 'strong20x1020240327_024931_sdl_s1_sg1';
  inFilePrePre = 'string20x1020240330_142942_sdl_s8_sg2';
  inFilePrePre = 'string20x1020240330_173800_sdl_s9_sg2';
  inFilePrePre = 'string20x1020240330_174216_sdl_s10_sg1';
  inFilePre = [inFilePrePre]; 
  
  inPath = fullfile(vtData, 'brachial_vein');
  inFilePrePre = 'brachial_vein20240405_155355_sdl_s12_sg1';
  inFilePre = [inFilePrePre]; 
  
  % sg1-4 have signal
  inPath = fullfile(vtData, 'Apr-11-pig');
  inFilePrePre = '120240411_215120_sdl_s1_sg';
  sgSel = [1:4]; 
  % x,z ROI that encompasses whole vessel segment of interest
  xzRect = [-2.1305    7.6530
            -1.5098    9.1309];
end


if 0
 inPath = fullfile(vtData, '20240501_Survival');
 % 03_withnoises20240501_130351_sdl_s3_sg bad pulse noise
 inFilePrePre = '0420240501_133240_sdl_s1_sg'; % quite bad
 inFilePrePre = '0520240501_133445_sdl_s2_sg'; % bad
 inPath = fullfile(vtData, '20240501_Survival_post');
 inFilePrePre = 'distal_injury120240501_152544_sdl_s8_sg'; % nothing
 inFilePrePre = 'distal_injury220240501_152637_sdl_s9_sg'; % nothing
 inFilePrePre = 'healthyref_0120240501_144212_sdl_s3_sg'; % bad
 inFilePrePre = 'healthyref_0220240501_144359_sdl_s4_sg';
 inFilePrePre = 'mid_injury120240501_151414_sdl_s5_sg'; % some Doppler, but weak
 inFilePrePre = 'mid_injury220240501_151510_sdl_s6_sg'; % some Doppler, but weak
 inFilePrePre = 'mid_injury320240501_151616_sdl_s7_sg'; % some, very weak
 sgSel = [1];
 xzRect = [-2.1305    7.6530
           -1.5098    9.1309];
end

if 0
  inPath = fullfile(vtData, '0627MUSICPIG');
  inFilePrePre = 'FUS_pig_afterkataminepreFUS20240627_194601_sdl_s1_sg';
  inFilePre = [inFilePrePre]; 
  sgSel = [1:6];
  %sgSel = [1];
  % restrict search for max flow pixels to area with vessel
  xzRect = [0 6
            1.3 8.21];
  
end

if 0
  inPath = fullfile(vtData, '06272024MUSIC');
  inFilePrePre = 'PostPostFUS20240627_200144_sdl_s1_sg';
  inFilePre = [inFilePrePre]; 
  sgSel = [1];
  %sgSel = [1];
  % restrict search for max flow pixels to area with vessel
  xzRect = [0 6
            1.3 8.21];
  
end


if 0
  inPath = fullfile(vtData, '20240718');
  inFilePrePre = 'PreFUS20240718_165558_sdl_s1_sg';
  inFilePre = [inFilePrePre]; 
  sgSel = [1]; % not clear if signals are from artery or vein
               %  sgSel = [2]; % not clear if signals are from artery or vein
               %  sgSel = [3]; % not clear if signals are from artery or vein
               %  sgSel = [4]; % not clear if signals are from artery or vein
               %  sgSel = [5]; % not clear if signals are from artery or vein
               %  sgSel = [6]; % not clear if signals are from artery or vein
  sgSel = [1:6]; % not clear if signals are from artery or vein          
  sgSel = [1];
  % restrict search for max flow pixels to area with vessel

  % for including max power cluster that may not be sulcal
  if 1
    xzRect = [0.75   6
              0.95   9.5]; % avoid outlier neat x=0, z=8, or that could be the real artery site.
                     % not clear
                     %  xzRect = [-2  6
                     %2   11]; % avoid outlier neat x=0, z=8, or that could be the real artery site.
                     % not clear
  else
      xzRect = [0 6
                0.25 11];

      %      [17,2]
  end
  
  
end


if 0
  inPath = fullfile(vtData, '20240718');
  inFilePrePre = 'PostFUS20240718_170124_sdl_s2_sg';
  inFilePre = [inFilePrePre]; 
  sgSel = [1:9]; % not clear if signals are from artery or vein
               % sgSel = [2]; % not clear if signals are from artery or vein
               %  sgSel = [3]; % not clear if signals are from artery or vein
               %  sgSel = [4]; % not clear if signals are from artery or vein
               %  sgSel = [5]; % not clear if signals are from artery or vein
               %  sgSel = [6]; % not clear if signals are from artery or vein
               %  sgSel = [1:6]; % not clear if signals are from artery or vein          
                 sgSel = [1];
  % restrict search for max flow pixels to area with vessel
  xzRect = [0.66 7
            2   11]; % avoid outlier neat x=0, z=8, or that could be the real artery site.
                     % not clear
                     % signals 
  xzRect = [0.9   8
            1.0  9.0];
end


if 0
  inPath = fullfile(vtData, '20240729');
  inFilePrePre = 'FUS1pre20240729_162348_sdl_s1_sg'; % nice, clear
  inFilePrePre = 'FUS1pre(opt)20240729_162443_sdl_s2_sg'; % nice, clear  
  inFilePre = [inFilePrePre]; 
  sgSel = [1:27]; 
  %sgSel = [1];
  % restrict search for max flow pixels to area with vessel
  xzRect = [-6   4
            6   11];
end


if 0
  inPath = fullfile(vtData, '20240729');
  inFilePrePre = 'FUS1pre(opt)20240729_162443_sdl_s2_sg'; % nice, clear  
  inFilePre = [inFilePrePre]; 
  sgSel = [1:27]; 
  %sgSel = [1];
  % restrict search for max flow pixels to area with vessel
  xzRect = [-6 4
            6  8]; % remove lower edge region from consideration
end


if 0
  inPath = fullfile(vtData, '20240729');
  inFilePrePre = 'FUS1post20240729_163730_sdl_s3_sg'; % good
  inFilePre = [inFilePrePre]; 
  sgSel = [1:15]; 
  %sgSel = [1];
  % restrict search for max flow pixels to area with vessel
  xzRect = [-6 5.4
            6  8]; % remove lower edge region from consideration
end


if 0
  inPath = fullfile(vtData, '20240729');
  inFilePrePre = 'FUS2post20240729_164533_sdl_s4_sg'; % good
  inFilePre = [inFilePrePre]; 
  sgSel = [1:15]; 
  %sgSel = [1];
  % restrict search for max flow pixels to area with vessel
  xzRect = [-6 5.4
            6  8]; % remove lower edge region from consideration
end


if 0
  inPath = fullfile(vtData, '20240729');
  inFilePrePre = 'FUS3post20240729_165341_sdl_s5_sg'; % good
  inFilePre = [inFilePrePre]; 
  sgSel = [1:14]; 
  %sgSel = [1];
  % restrict search for max flow pixels to area with vessel
  xzRect = [-6 5.4
            6  7.5]; % remove lower edge region from consideration
end


if 0
  inPath = fullfile(vtData, '20240729');
  inFilePrePre = 'FUS4post20240729_170139_sdl_s6_sg'; % good
  inFilePre = [inFilePrePre]; 
  sgSel = [1:26]; 
  %sgSel = [1];
  % restrict search for max flow pixels to area with vessel
  xzRect = [-6 5.4
            6  7.5]; % remove lower edge region from consideration
end

if 1
  inPath = fullfile(vtData, '0815Survival')
  inFilePrePre = 'implant_injurtsite120240815_132555_sdl_s1_sg'; sgSel=1;
  inFilePrePre = 'implant_injurtsite120240815_135106_sdl_s2_sg'; sgSel=1;
  inFilePrePre = 'healthy_pre20240815_140736_sdl_s3_sg'; sgSel=1:3;
  inFilePre = [inFilePrePre]; 
  %sgSel = [1];
  % restrict search for max flow pixels to area with vessel
  xzRect = [-6 5.4
            6  7.5]; % remove lower edge region from consideration
end


if 0
  inPath = fullfile(vtData, '815Survival')
  inFilePrePre = 'Post _injury_injsite20240815_164001_sdl_s1_sg'; sgSel=1;
  inFilePrePre = 'Post _injury_injsite20240815_164059_sdl_s2_sg'; sgSel=1;
  inFilePrePre = 'Post _injury_healthy20240815_164221_sdl_s3_sg'; sgSel=1:2;
  inFilePrePre = 'Post _injury_injsite20240815_164442_sdl_s4_sg'; sgSel=1:2;    
  inFilePrePre = 'Post _injury_injsite20240815_164524_sdl_s5_sg'; sgSel=1:2;
  
  
  inFilePre = [inFilePrePre]; 
  %sgSel = [1];
  % restrict search for max flow pixels to area with vessel
  xzRect = [-6 5.4
            6  7.5]; % remove lower edge region from consideration
end


if 0
  inPath = fullfile(vtData, '816dayone_survival')
  inFilePrePre = 'Healthy_array20240816_101001_sdl_s1_sg'; sgSel=1:2;
  inFilePrePre = 'Healthy_array20240816_101001_sdl_s1_sg'; sgSel=1:2;

  
  inFilePrePre = 'Healthy_array20240816_101104_sdl_s2_sg'; sgSel=1:4;
  inFilePrePre = 'Healthy_array20240816_101352_sdl_s3_sg'; sgSel=1;  
  
  inFilePrePre = 'healhty_retry120240816_105009_sdl_s13_sg'; sgSel=1:2;  
  inFilePrePre = 'healhty_retry120240816_105009_sdl_s13_sg'; sgSel=1:2;

  
  
  inFilePre = [inFilePrePre]; 
  %sgSel = [1];
  % restrict search for max flow pixels to area with vessel
  xzRect = [-6 5.4
            6  7.5]; % remove lower edge region from consideration
end


 % r=1:10; c=1:5;
 %r=10:11; c=5;
if 0
  if 0
  x_mm = 0.5;
  z_mm = 48;  
  xWid_mm=2;
  zWid_mm=2;
  xzOverride=0;
  imPlotInd = 200;
  end
  
end
