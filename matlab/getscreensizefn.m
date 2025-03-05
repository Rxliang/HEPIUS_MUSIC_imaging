function screensize = getscreensizefn

ge = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment;
gd = ge.getDefaultScreenDevice;
screensize = [gd.getDisplayMode.getWidth gd.getDisplayMode.getHeight];