vtData = getenv('MUSIC_DATA');
music_home = getenv('MUSIC_HOME');

in_file1 = fullfile(vtData, '816dayone_survival', 'iq2sono_fft_Healthy_array120240816_101104_sdl_s2_1_2_3_4_ind380.mat');

in_file2 = fullfile(vtData, '816dayone_survival', 'iq2sono_fft_Healthy_array120240816_101104_sdl_s2_1_2_3_4_ind460.mat');


s1 = load(in_file1);
s2 = load(in_file2);

figure(ff)
clf
  %  imagesc(timeAxFull_s, sonoSession.fAx,  log10(sonoImFull), [10 20]);
  %  im_lim = [13 16];
  im_lim = [];
  imagesc(timeAxSel_s, sonoSession.fAx,  log10(s1.sonoImSel) - log10(s2.sonoImSel))%, im_lim);  
  colormap(gray)
  axis xy
  xlabel('time (s)')
  ylabel('frequency (Hz)')
  ind_rc = ij2ind(r,c, sz(1));
  title(['Gate at row ' num2str(r) ', column ' num2str(c) ', ind ' num2str(ind_rc)]);
  % title(tldc, replacechar(inFilePrePre, '_', '\_'))




