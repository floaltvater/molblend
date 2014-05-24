import bpy
import os

wm = bpy.context.window_manager
kc = wm.keyconfigs.new("MolBlend")

# Map 3D View
km = kc.keymaps.new('3D View', space_type='VIEW_3D', region_type='WINDOW', modal=False)

try:
    kmitems = km.keymap_items
    kmi = kmitems.new('mb.left_click', 'LEFTMOUSE', 'PRESS')
    kmi.properties.shift = False
    kmi.properties.ctrl = False
    kmi.properties.alt = False
    kmi = kmitems.new('mb.left_click', 'LEFTMOUSE', 'PRESS', shift=True)
    kmi.properties.shift = True
    kmi.properties.ctrl = False
    kmi.properties.alt = False
    kmi = kmitems.new('mb.left_click', 'LEFTMOUSE', 'PRESS', ctrl=True)
    kmi.properties.shift = False
    kmi.properties.ctrl = True
    kmi.properties.alt = False
    kmi = kmitems.new('mb.left_click', 'LEFTMOUSE', 'PRESS', alt=True)
    kmi.properties.shift = False
    kmi.properties.ctrl = False
    kmi.properties.alt = True
    kmi = kmitems.new('mb.left_click', 'LEFTMOUSE', 'PRESS', shift=True, ctrl=True)
    kmi.properties.shift = True
    kmi.properties.ctrl = True
    kmi.properties.alt = False
    kmi = kmitems.new('mb.left_click', 'LEFTMOUSE', 'PRESS', shift=True, alt=True)
    kmi.properties.shift = True
    kmi.properties.ctrl = False
    kmi.properties.alt = True
    kmi = kmitems.new('mb.left_click', 'LEFTMOUSE', 'PRESS', ctrl=True, alt=True)
    kmi.properties.shift = False
    kmi.properties.ctrl = True
    kmi.properties.alt = True
    kmi = kmitems.new('mb.left_click', 'LEFTMOUSE', 'PRESS', shift=True, ctrl=True, alt=True)
    kmi.properties.shift = True
    kmi.properties.ctrl = True
    kmi.properties.alt = True
    kmi = kmitems.new('mb.right_click', 'RIGHTMOUSE', 'PRESS')
    #kmi.properties.shift = False
    #kmi.properties.ctrl = False
    #kmi.properties.alt = False
    #context.window_manager.mb.added_kmitems.add(kmi.id)

    ### DEFAULT keymap

    #kmi = km.keymap_items.new('view3d.manipulator', 'LEFTMOUSE', 'PRESS', any=True)
    #kmi.properties.release_confirm = True
    #kmi = km.keymap_items.new('view3d.cursor3d', 'ACTIONMOUSE', 'PRESS')
    kmi = km.keymap_items.new('view3d.rotate', 'MIDDLEMOUSE', 'PRESS')
    kmi = km.keymap_items.new('view3d.move', 'MIDDLEMOUSE', 'PRESS', shift=True)
    kmi = km.keymap_items.new('view3d.zoom', 'MIDDLEMOUSE', 'PRESS', ctrl=True)
    kmi = km.keymap_items.new('view3d.dolly', 'MIDDLEMOUSE', 'PRESS', shift=True, ctrl=True)
    kmi = km.keymap_items.new('view3d.view_selected', 'NUMPAD_PERIOD', 'PRESS', ctrl=True)
    kmi.properties.use_all_regions = True
    kmi = km.keymap_items.new('view3d.view_selected', 'NUMPAD_PERIOD', 'PRESS')
    kmi.properties.use_all_regions = False
    kmi = km.keymap_items.new('view3d.view_lock_to_active', 'NUMPAD_PERIOD', 'PRESS', shift=True)
    kmi = km.keymap_items.new('view3d.view_lock_clear', 'NUMPAD_PERIOD', 'PRESS', alt=True)
    kmi = km.keymap_items.new('view3d.fly', 'F', 'PRESS', shift=True)
    kmi = km.keymap_items.new('view3d.smoothview', 'TIMER1', 'ANY', any=True)
    kmi = km.keymap_items.new('view3d.rotate', 'TRACKPADPAN', 'ANY')
    kmi = km.keymap_items.new('view3d.rotate', 'MOUSEROTATE', 'ANY')
    kmi = km.keymap_items.new('view3d.move', 'TRACKPADPAN', 'ANY', shift=True)
    kmi = km.keymap_items.new('view3d.zoom', 'TRACKPADZOOM', 'ANY')
    kmi = km.keymap_items.new('view3d.zoom', 'TRACKPADPAN', 'ANY', ctrl=True)
    kmi = km.keymap_items.new('view3d.zoom', 'NUMPAD_PLUS', 'PRESS')
    kmi.properties.delta = 1
    kmi = km.keymap_items.new('view3d.zoom', 'NUMPAD_MINUS', 'PRESS')
    kmi.properties.delta = -1
    kmi = km.keymap_items.new('view3d.zoom', 'EQUAL', 'PRESS', ctrl=True)
    kmi.properties.delta = 1
    kmi = km.keymap_items.new('view3d.zoom', 'MINUS', 'PRESS', ctrl=True)
    kmi.properties.delta = -1
    kmi = km.keymap_items.new('view3d.zoom', 'WHEELINMOUSE', 'PRESS')
    kmi.properties.delta = 1
    kmi = km.keymap_items.new('view3d.zoom', 'WHEELOUTMOUSE', 'PRESS')
    kmi.properties.delta = -1
    kmi = km.keymap_items.new('view3d.dolly', 'NUMPAD_PLUS', 'PRESS', shift=True)
    kmi.properties.delta = 1
    kmi = km.keymap_items.new('view3d.dolly', 'NUMPAD_MINUS', 'PRESS', shift=True)
    kmi.properties.delta = -1
    kmi = km.keymap_items.new('view3d.dolly', 'EQUAL', 'PRESS', shift=True, ctrl=True)
    kmi.properties.delta = 1
    kmi = km.keymap_items.new('view3d.dolly', 'MINUS', 'PRESS', shift=True, ctrl=True)
    kmi.properties.delta = -1
    kmi = km.keymap_items.new('view3d.zoom_camera_1_to_1', 'NUMPAD_ENTER', 'PRESS', shift=True)
    kmi = km.keymap_items.new('view3d.view_center_camera', 'HOME', 'PRESS')
    kmi = km.keymap_items.new('view3d.view_center_lock', 'HOME', 'PRESS')
    kmi = km.keymap_items.new('view3d.view_center_cursor', 'HOME', 'PRESS', alt=True)
    kmi = km.keymap_items.new('view3d.view_center_pick', 'F', 'PRESS', alt=True)
    kmi = km.keymap_items.new('view3d.view_all', 'HOME', 'PRESS')
    kmi.properties.center = False
    kmi = km.keymap_items.new('view3d.view_all', 'HOME', 'PRESS', ctrl=True)
    kmi.properties.use_all_regions = True
    kmi.properties.center = False
    kmi = km.keymap_items.new('view3d.view_all', 'C', 'PRESS', shift=True)
    kmi.properties.center = True
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NUMPAD_0', 'PRESS')
    kmi.properties.type = 'CAMERA'
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NUMPAD_1', 'PRESS')
    kmi.properties.type = 'FRONT'
    kmi = km.keymap_items.new('view3d.view_orbit', 'NUMPAD_2', 'PRESS')
    kmi.properties.type = 'ORBITDOWN'
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NUMPAD_3', 'PRESS')
    kmi.properties.type = 'RIGHT'
    kmi = km.keymap_items.new('view3d.view_orbit', 'NUMPAD_4', 'PRESS')
    kmi.properties.type = 'ORBITLEFT'
    kmi = km.keymap_items.new('view3d.view_persportho', 'NUMPAD_5', 'PRESS')
    kmi = km.keymap_items.new('view3d.view_orbit', 'NUMPAD_6', 'PRESS')
    kmi.properties.type = 'ORBITRIGHT'
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NUMPAD_7', 'PRESS')
    kmi.properties.type = 'TOP'
    kmi = km.keymap_items.new('view3d.view_orbit', 'NUMPAD_8', 'PRESS')
    kmi.properties.type = 'ORBITUP'
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NUMPAD_1', 'PRESS', ctrl=True)
    kmi.properties.type = 'BACK'
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NUMPAD_3', 'PRESS', ctrl=True)
    kmi.properties.type = 'LEFT'
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NUMPAD_7', 'PRESS', ctrl=True)
    kmi.properties.type = 'BOTTOM'
    kmi = km.keymap_items.new('view3d.view_pan', 'NUMPAD_2', 'PRESS', ctrl=True)
    kmi.properties.type = 'PANDOWN'
    kmi = km.keymap_items.new('view3d.view_pan', 'NUMPAD_4', 'PRESS', ctrl=True)
    kmi.properties.type = 'PANLEFT'
    kmi = km.keymap_items.new('view3d.view_pan', 'NUMPAD_6', 'PRESS', ctrl=True)
    kmi.properties.type = 'PANRIGHT'
    kmi = km.keymap_items.new('view3d.view_pan', 'NUMPAD_8', 'PRESS', ctrl=True)
    kmi.properties.type = 'PANUP'
    kmi = km.keymap_items.new('view3d.view_roll', 'NUMPAD_4', 'PRESS', shift=True)
    kmi.properties.angle = -0.2617993950843811
    kmi = km.keymap_items.new('view3d.view_roll', 'NUMPAD_6', 'PRESS', shift=True)
    kmi.properties.angle = 0.2617993950843811
    kmi = km.keymap_items.new('view3d.view_pan', 'WHEELUPMOUSE', 'PRESS', ctrl=True)
    kmi.properties.type = 'PANRIGHT'
    kmi = km.keymap_items.new('view3d.view_pan', 'WHEELDOWNMOUSE', 'PRESS', ctrl=True)
    kmi.properties.type = 'PANLEFT'
    kmi = km.keymap_items.new('view3d.view_pan', 'WHEELUPMOUSE', 'PRESS', shift=True)
    kmi.properties.type = 'PANUP'
    kmi = km.keymap_items.new('view3d.view_pan', 'WHEELDOWNMOUSE', 'PRESS', shift=True)
    kmi.properties.type = 'PANDOWN'
    kmi = km.keymap_items.new('view3d.view_orbit', 'WHEELUPMOUSE', 'PRESS', ctrl=True, alt=True)
    kmi.properties.type = 'ORBITLEFT'
    kmi = km.keymap_items.new('view3d.view_orbit', 'WHEELDOWNMOUSE', 'PRESS', ctrl=True, alt=True)
    kmi.properties.type = 'ORBITRIGHT'
    kmi = km.keymap_items.new('view3d.view_orbit', 'WHEELUPMOUSE', 'PRESS', shift=True, alt=True)
    kmi.properties.type = 'ORBITUP'
    kmi = km.keymap_items.new('view3d.view_orbit', 'WHEELDOWNMOUSE', 'PRESS', shift=True, alt=True)
    kmi.properties.type = 'ORBITDOWN'
    kmi = km.keymap_items.new('view3d.view_roll', 'WHEELUPMOUSE', 'PRESS', shift=True, ctrl=True)
    kmi.properties.angle = -0.2617993950843811
    kmi = km.keymap_items.new('view3d.view_roll', 'WHEELDOWNMOUSE', 'PRESS', shift=True, ctrl=True)
    kmi.properties.angle = 0.2617993950843811
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NUMPAD_1', 'PRESS', shift=True)
    kmi.properties.type = 'FRONT'
    kmi.properties.align_active = True
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NUMPAD_3', 'PRESS', shift=True)
    kmi.properties.type = 'RIGHT'
    kmi.properties.align_active = True
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NUMPAD_7', 'PRESS', shift=True)
    kmi.properties.type = 'TOP'
    kmi.properties.align_active = True
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NUMPAD_1', 'PRESS', shift=True, ctrl=True)
    kmi.properties.type = 'BACK'
    kmi.properties.align_active = True
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NUMPAD_3', 'PRESS', shift=True, ctrl=True)
    kmi.properties.type = 'LEFT'
    kmi.properties.align_active = True
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NUMPAD_7', 'PRESS', shift=True, ctrl=True)
    kmi.properties.type = 'BOTTOM'
    kmi.properties.align_active = True
    kmi = km.keymap_items.new('view3d.localview', 'NUMPAD_SLASH', 'PRESS')
    kmi = km.keymap_items.new('view3d.ndof_orbit_zoom', 'NDOF_MOTION', 'ANY')
    kmi = km.keymap_items.new('view3d.ndof_orbit', 'NDOF_MOTION', 'ANY', ctrl=True)
    kmi = km.keymap_items.new('view3d.ndof_pan', 'NDOF_MOTION', 'ANY', shift=True)
    kmi = km.keymap_items.new('view3d.ndof_all', 'NDOF_MOTION', 'ANY', shift=True, ctrl=True)
    kmi = km.keymap_items.new('view3d.view_selected', 'NDOF_BUTTON_FIT', 'PRESS')
    kmi.properties.use_all_regions = False
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NDOF_BUTTON_FRONT', 'PRESS')
    kmi.properties.type = 'FRONT'
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NDOF_BUTTON_BACK', 'PRESS')
    kmi.properties.type = 'BACK'
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NDOF_BUTTON_LEFT', 'PRESS')
    kmi.properties.type = 'LEFT'
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NDOF_BUTTON_RIGHT', 'PRESS')
    kmi.properties.type = 'RIGHT'
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NDOF_BUTTON_TOP', 'PRESS')
    kmi.properties.type = 'TOP'
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NDOF_BUTTON_BOTTOM', 'PRESS')
    kmi.properties.type = 'BOTTOM'
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NDOF_BUTTON_FRONT', 'PRESS', shift=True)
    kmi.properties.type = 'FRONT'
    kmi.properties.align_active = True
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NDOF_BUTTON_RIGHT', 'PRESS', shift=True)
    kmi.properties.type = 'RIGHT'
    kmi.properties.align_active = True
    kmi = km.keymap_items.new('view3d.viewnumpad', 'NDOF_BUTTON_TOP', 'PRESS', shift=True)
    kmi.properties.type = 'TOP'
    kmi.properties.align_active = True
    kmi = km.keymap_items.new('view3d.layers', 'ACCENT_GRAVE', 'PRESS')
    kmi.properties.nr = 0
    kmi = km.keymap_items.new('view3d.layers', 'ONE', 'PRESS', any=True)
    kmi.properties.nr = 1
    kmi = km.keymap_items.new('view3d.layers', 'TWO', 'PRESS', any=True)
    kmi.properties.nr = 2
    kmi = km.keymap_items.new('view3d.layers', 'THREE', 'PRESS', any=True)
    kmi.properties.nr = 3
    kmi = km.keymap_items.new('view3d.layers', 'FOUR', 'PRESS', any=True)
    kmi.properties.nr = 4
    kmi = km.keymap_items.new('view3d.layers', 'FIVE', 'PRESS', any=True)
    kmi.properties.nr = 5
    kmi = km.keymap_items.new('view3d.layers', 'SIX', 'PRESS', any=True)
    kmi.properties.nr = 6
    kmi = km.keymap_items.new('view3d.layers', 'SEVEN', 'PRESS', any=True)
    kmi.properties.nr = 7
    kmi = km.keymap_items.new('view3d.layers', 'EIGHT', 'PRESS', any=True)
    kmi.properties.nr = 8
    kmi = km.keymap_items.new('view3d.layers', 'NINE', 'PRESS', any=True)
    kmi.properties.nr = 9
    kmi = km.keymap_items.new('view3d.layers', 'ZERO', 'PRESS', any=True)
    kmi.properties.nr = 10
    kmi = km.keymap_items.new('wm.context_toggle_enum', 'Z', 'PRESS')
    kmi.properties.data_path = 'space_data.viewport_shade'
    kmi.properties.value_1 = 'SOLID'
    kmi.properties.value_2 = 'WIREFRAME'
    kmi = km.keymap_items.new('wm.context_toggle_enum', 'Z', 'PRESS', alt=True)
    kmi.properties.data_path = 'space_data.viewport_shade'
    kmi.properties.value_1 = 'SOLID'
    kmi.properties.value_2 = 'TEXTURED'
    kmi = km.keymap_items.new('wm.context_toggle_enum', 'Z', 'PRESS', shift=True)
    kmi.properties.data_path = 'space_data.viewport_shade'
    kmi.properties.value_1 = 'SOLID'
    kmi.properties.value_2 = 'RENDERED'
    #kmi = km.keymap_items.new('view3d.select', 'SELECTMOUSE', 'PRESS')
    #kmi.properties.extend = False
    #kmi.properties.deselect = False
    #kmi.properties.toggle = False
    #kmi.properties.center = False
    #kmi.properties.enumerate = False
    #kmi.properties.object = False
    #kmi = km.keymap_items.new('view3d.select', 'SELECTMOUSE', 'PRESS', shift=True)
    #kmi.properties.extend = False
    #kmi.properties.deselect = False
    #kmi.properties.toggle = True
    #kmi.properties.center = False
    #kmi.properties.enumerate = False
    #kmi.properties.object = False
    #kmi = km.keymap_items.new('view3d.select', 'SELECTMOUSE', 'PRESS', ctrl=True)
    #kmi.properties.extend = False
    #kmi.properties.deselect = False
    #kmi.properties.toggle = False
    #kmi.properties.center = True
    #kmi.properties.enumerate = False
    #kmi.properties.object = True
    #kmi = km.keymap_items.new('view3d.select', 'SELECTMOUSE', 'PRESS', alt=True)
    #kmi.properties.extend = False
    #kmi.properties.deselect = False
    #kmi.properties.toggle = False
    #kmi.properties.center = False
    #kmi.properties.enumerate = True
    #kmi.properties.object = False
    #kmi = km.keymap_items.new('view3d.select', 'SELECTMOUSE', 'PRESS', shift=True, ctrl=True)
    #kmi.properties.extend = True
    #kmi.properties.deselect = False
    #kmi.properties.toggle = True
    #kmi.properties.center = True
    #kmi.properties.enumerate = False
    #kmi.properties.object = False
    #kmi = km.keymap_items.new('view3d.select', 'SELECTMOUSE', 'PRESS', ctrl=True, alt=True)
    #kmi.properties.extend = False
    #kmi.properties.deselect = False
    #kmi.properties.toggle = False
    #kmi.properties.center = True
    #kmi.properties.enumerate = True
    #kmi.properties.object = False
    #kmi = km.keymap_items.new('view3d.select', 'SELECTMOUSE', 'PRESS', shift=True, alt=True)
    #kmi.properties.extend = False
    #kmi.properties.deselect = False
    #kmi.properties.toggle = True
    #kmi.properties.center = False
    #kmi.properties.enumerate = True
    #kmi.properties.object = False
    #kmi = km.keymap_items.new('view3d.select', 'SELECTMOUSE', 'PRESS', shift=True, ctrl=True, alt=True)
    #kmi.properties.extend = False
    #kmi.properties.deselect = False
    #kmi.properties.toggle = True
    #kmi.properties.center = True
    #kmi.properties.enumerate = True
    #kmi.properties.object = False
    kmi = km.keymap_items.new('view3d.select_border', 'B', 'PRESS')
    kmi = km.keymap_items.new('view3d.select_lasso', 'EVT_TWEAK_A', 'ANY', ctrl=True)
    kmi.properties.deselect = False
    kmi = km.keymap_items.new('view3d.select_lasso', 'EVT_TWEAK_A', 'ANY', shift=True, ctrl=True)
    kmi.properties.deselect = True
    kmi = km.keymap_items.new('view3d.select_circle', 'C', 'PRESS')
    kmi = km.keymap_items.new('view3d.clip_border', 'B', 'PRESS', alt=True)
    kmi = km.keymap_items.new('view3d.zoom_border', 'B', 'PRESS', shift=True)
    kmi = km.keymap_items.new('view3d.render_border', 'B', 'PRESS', shift=True)
    kmi.properties.camera_only = True
    kmi = km.keymap_items.new('view3d.render_border', 'B', 'PRESS', ctrl=True)
    kmi.properties.camera_only = False
    kmi = km.keymap_items.new('view3d.clear_render_border', 'B', 'PRESS', ctrl=True, alt=True)
    kmi = km.keymap_items.new('view3d.camera_to_view', 'NUMPAD_0', 'PRESS', ctrl=True, alt=True)
    kmi = km.keymap_items.new('view3d.object_as_camera', 'NUMPAD_0', 'PRESS', ctrl=True)
    kmi = km.keymap_items.new('wm.call_menu', 'S', 'PRESS', shift=True)
    kmi.properties.name = 'VIEW3D_MT_snap'
    kmi = km.keymap_items.new('view3d.copybuffer', 'C', 'PRESS', ctrl=True)
    kmi = km.keymap_items.new('view3d.pastebuffer', 'V', 'PRESS', ctrl=True)
    kmi = km.keymap_items.new('wm.context_set_enum', 'COMMA', 'PRESS')
    kmi.properties.data_path = 'space_data.pivot_point'
    kmi.properties.value = 'BOUNDING_BOX_CENTER'
    kmi = km.keymap_items.new('wm.context_set_enum', 'COMMA', 'PRESS', ctrl=True)
    kmi.properties.data_path = 'space_data.pivot_point'
    kmi.properties.value = 'MEDIAN_POINT'
    kmi = km.keymap_items.new('wm.context_toggle', 'COMMA', 'PRESS', alt=True)
    kmi.properties.data_path = 'space_data.use_pivot_point_align'
    kmi = km.keymap_items.new('wm.context_toggle', 'SPACE', 'PRESS', ctrl=True)
    kmi.properties.data_path = 'space_data.show_manipulator'
    kmi = km.keymap_items.new('wm.context_set_enum', 'PERIOD', 'PRESS')
    kmi.properties.data_path = 'space_data.pivot_point'
    kmi.properties.value = 'CURSOR'
    kmi = km.keymap_items.new('wm.context_set_enum', 'PERIOD', 'PRESS', ctrl=True)
    kmi.properties.data_path = 'space_data.pivot_point'
    kmi.properties.value = 'INDIVIDUAL_ORIGINS'
    kmi = km.keymap_items.new('wm.context_set_enum', 'PERIOD', 'PRESS', alt=True)
    kmi.properties.data_path = 'space_data.pivot_point'
    kmi.properties.value = 'ACTIVE_ELEMENT'
    kmi = km.keymap_items.new('transform.translate', 'G', 'PRESS')
    kmi = km.keymap_items.new('transform.translate', 'EVT_TWEAK_S', 'ANY')
    kmi = km.keymap_items.new('transform.rotate', 'R', 'PRESS')
    kmi = km.keymap_items.new('transform.resize', 'S', 'PRESS')
    kmi = km.keymap_items.new('transform.warp', 'W', 'PRESS', shift=True)
    kmi = km.keymap_items.new('transform.tosphere', 'S', 'PRESS', shift=True, alt=True)
    kmi = km.keymap_items.new('transform.shear', 'S', 'PRESS', shift=True, ctrl=True, alt=True)
    kmi = km.keymap_items.new('transform.select_orientation', 'SPACE', 'PRESS', alt=True)
    kmi = km.keymap_items.new('transform.create_orientation', 'SPACE', 'PRESS', ctrl=True, alt=True)
    kmi.properties.use = True
    kmi = km.keymap_items.new('transform.mirror', 'M', 'PRESS', ctrl=True)
    kmi = km.keymap_items.new('wm.context_toggle', 'TAB', 'PRESS', shift=True)
    kmi.properties.data_path = 'tool_settings.use_snap'
    kmi = km.keymap_items.new('wm.context_menu_enum', 'TAB', 'PRESS', shift=True, ctrl=True)
    kmi.properties.data_path = 'tool_settings.snap_element'
    kmi = km.keymap_items.new('transform.translate', 'T', 'PRESS', shift=True)
    kmi.properties.texture_space = True
    kmi = km.keymap_items.new('transform.resize', 'T', 'PRESS', shift=True, alt=True)
    kmi.properties.texture_space = True
    kmi = km.keymap_items.new('transform.skin_resize', 'A', 'PRESS', ctrl=True)
except AttributeError:
    print(kmi.name)
