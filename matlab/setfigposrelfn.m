function pos_abs = setfigposrelfn(gcf, pos_rel)

screen_size = get(0,'ScreenSize');
pos_abs =  pos_rel.*screen_size;
set(gcf, 'position', pos_abs);

end

