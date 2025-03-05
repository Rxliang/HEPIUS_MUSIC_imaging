from dearpygui import dearpygui


dearpygui.create_context()
dearpygui.setup_dearpygui()
dearpygui.create_viewport()
dearpygui.show_viewport()


#font_fpath = r'Inter-Regular.ttf'   # | update with local font file!
#font_name  = 'Inter-Regular'        # |

font_fpath = r'c:/windows/fonts/calibri.ttf'   # | update with local font file!
font_name  = 'Calibri-Regular'        # |


with dearpygui.window(): ... # needed to populate the item registry

with dearpygui.font_registry(label="Oldest") as fr_oldest:
    for i in range(10, 20, 2):
        dearpygui.add_font(font_fpath, i, label=f'{font_name} [{i:.1f}]')

with dearpygui.font_registry(label="Newest") as fr_newest:
    for i in range(20, 30, 2):
        dearpygui.add_font(font_fpath, i, label=f'{font_name} [{i:.1f}]')

# will still reflect the built-in ProggyClean 13px font
dearpygui.bind_font(dearpygui.last_item())


dearpygui.show_font_manager()
dearpygui.show_item_registry()
dearpygui.start_dearpygui()
