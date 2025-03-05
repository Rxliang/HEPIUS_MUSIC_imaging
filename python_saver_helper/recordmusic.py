#!C:\Program Files\Python312\pythonw.exe
import dearpygui.dearpygui as dpg
import pyautogui
import math
from setsemaon import setsemaon
from setsemaoff import setsemaoff

step_no = 0;

class LayoutHelper:

    def __init__(self):
        self.table_id = dpg.add_table(header_row=False, policy=dpg.mvTable_SizingStretchProp)
        self.stage_id = dpg.add_stage()
        dpg.push_container_stack(self.stage_id)
        
    def add_widget(self, uuid, percentage):
        dpg.add_table_column(init_width_or_weight=percentage/100.0, parent=self.table_id)
        dpg.set_item_width(uuid, -1)
        return uuid
    
    def submit(self):
        dpg.pop_container_stack() # pop stage
        with dpg.table_row(parent=self.table_id):
            dpg.unstage(self.stage_id)
            
BACKGROUND_COLOUR = (0,0,0)

sz = pyautogui.size()
print(str(sz))
ui_wid = math.ceil(sz[0]/6)
ui_ht = math.ceil(sz[1]/3)

def button_callback(sender, app_data, user_data):
  global step_no
  # Unpack the user_data that is currently associated with the button
  state, enabled_theme, disabled_theme, subject, session, filetxt = user_data
  # Flip the state
  state = not state
  # Apply the appropriate theme
  dpg.bind_item_theme(sender, enabled_theme if state is True else disabled_theme)
  # Update the user_data associated with the button
  dpg.set_item_user_data(sender, (state, enabled_theme, disabled_theme, subject, session, filetxt))
  # this is the new state
  if state:
      setsemaoff()
      dpg.set_item_label(sender, "Start recording")
  else:
      subject_str = dpg.get_value(subject)      
      session_str = dpg.get_value(session)    
      step_no = step_no + 1
      filename = setsemaon(subject_str, session_str, step_no)
      print(filename)
      dpg.set_item_label(sender, "Stop recording")
      dpg.set_value(filetxt, 'output file: ' + filename)

  
dpg.create_context()

# add a font registry
with dpg.font_registry():
    # first argument ids the path to the .ttf or .otf file
    sz = pyautogui.size()
    fs = math.ceil(40*sz[0]/3840)
    default_font = dpg.add_font("c:/windows/fonts/calibri.ttf", fs)
  #  second_font = dpg.add_font("c:/windows/fonts/symbol.ttf", 5)
    second_font = dpg.add_font("c:/windows/fonts/arial.ttf", fs)
    
    
with dpg.theme() as global_theme:

    with dpg.theme_component(dpg.mvAll):
        dpg.add_theme_color(dpg.mvThemeCol_FrameBg, (70, 70, 70), category=dpg.mvThemeCat_Core)
        dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 5, category=dpg.mvThemeCat_Core)
        dpg.add_theme_color(dpg.mvThemeCol_Border, (1,0,0),  category=dpg.mvThemeCat_Core)
        dpg.add_theme_color(dpg.mvThemeCol_TitleBg, (1,0,0),  category=dpg.mvThemeCat_Core)
        dpg.add_theme_color(dpg.mvThemeCol_TitleBgActive, (1,0,0),  category=dpg.mvThemeCat_Core)
        dpg.add_theme_color(dpg.mvThemeCol_TitleBgCollapsed, (0.95,0,0.95),  category=dpg.mvThemeCat_Core)        

    with dpg.theme_component(dpg.mvInputInt):
        col = (0.1, 0.1, 0.1)
        dpg.add_theme_color(dpg.mvThemeCol_FrameBg, col, category=dpg.mvThemeCat_Core)
        dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 5, category=dpg.mvThemeCat_Core)


with dpg.theme() as enabled_theme:
    with dpg.theme_component(dpg.mvAll):
        dpg.add_theme_color(dpg.mvThemeCol_Text, (0, 255, 0),
                            category=dpg.mvThemeCat_Core)
# - This theme should make the label text on the button red

with dpg.theme() as disabled_theme:
    with dpg.theme_component(dpg.mvAll):
        dpg.add_theme_color(dpg.mvThemeCol_Text, (255, 0, 0),
                            category=dpg.mvThemeCat_Core)
        
        
dpg.bind_theme(global_theme)

dpg.create_viewport(title='MUSIC', width=ui_wid, height=ui_ht,
                    x_pos=sz[0]-ui_wid, y_pos=0, clear_color=(1,0,0))
#dpg.set_viewport_clear_color((0,0,0))
#dpg.set_viewport_title('fish')

dpg.setup_dearpygui()

with dpg.window(label="MUSIC recorder", width=ui_wid*0.96, height=ui_ht):
           
    with dpg.group(horizontal=True):
        dpg.add_text("session:") 
        s = dpg.add_input_text(tag="__session")

    dpg.add_text("")
    with dpg.group(horizontal=True):
        dpg.add_text("subject:") 
        sb = dpg.add_input_text(tag="__subject")
#        \u25a0
    dpg.add_text("")
    t = dpg.add_text("")
    dpg.add_text("")

    line1 = LayoutHelper()
    line1.add_widget(dpg.add_spacer(), 33.3)
    btn =line1.add_widget(
            dpg.add_button(tag='btn_tag',label="Start recording", callback=button_callback, user_data=(True, enabled_theme, disabled_theme, sb, s, t)), 60)                                 
    line1.add_widget(dpg.add_spacer(), 33.3)    
    line1.submit()

#    btn = dpg.get_value('btn_tag')
    
    #dpg.set_item_pos(btn, (dpg.get_viewport_width()/2, dpg.get_viewport_height()/2))
    dpg.bind_item_font(btn, second_font)
    dpg.bind_item_theme(btn, disabled_theme)

    
    # set font of specific widget
    dpg.bind_font(default_font)
#    dpg.bind_item_font(s, second_font)

    
#dpg.set_item_type_theme(dpg.mvButton, button_right)
dpg.show_viewport()
dpg.start_dearpygui()

dpg.destroy_context()
