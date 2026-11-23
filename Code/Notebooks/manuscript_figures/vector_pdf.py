"""Use installed librsvg/Cairo with a floating-point PDF viewport (no CLI rounding)."""
import ctypes as ct


def svg_to_pdf(svg_path, pdf_path, width_mm, height_mm):
    rsvg=ct.CDLL('librsvg-2.so.2')
    cairo=ct.CDLL('libcairo.so.2')
    gobject=ct.CDLL('libgobject-2.0.so.0')
    pointer=ct.c_void_p
    signatures={
        'cairo_pdf_surface_create':([ct.c_char_p,ct.c_double,ct.c_double],pointer),
        'cairo_create':([pointer],pointer),
        'cairo_surface_status':([pointer],ct.c_int),
        'cairo_destroy':([pointer],None),
        'cairo_surface_finish':([pointer],None),
        'cairo_surface_destroy':([pointer],None),
    }
    for name,(arguments,result) in signatures.items():
        function=getattr(cairo,name);function.argtypes=arguments;function.restype=result
    class Rectangle(ct.Structure):
        _fields_=[(name,ct.c_double) for name in ['x','y','width','height']]
    rsvg.rsvg_handle_new_from_file.argtypes=[ct.c_char_p,ct.POINTER(pointer)]
    rsvg.rsvg_handle_new_from_file.restype=pointer
    rsvg.rsvg_handle_set_dpi.argtypes=[pointer,ct.c_double]
    rsvg.rsvg_handle_render_document.argtypes=[pointer,pointer,ct.POINTER(Rectangle),ct.POINTER(pointer)]
    rsvg.rsvg_handle_render_document.restype=ct.c_int
    gobject.g_object_unref.argtypes=[pointer]
    handle=rsvg.rsvg_handle_new_from_file(str(svg_path).encode(),None)
    assert handle,svg_path
    width,height=width_mm*72/25.4,height_mm*72/25.4
    surface=cairo.cairo_pdf_surface_create(str(pdf_path).encode(),width,height)
    context=cairo.cairo_create(surface)
    try:
        rsvg.rsvg_handle_set_dpi(handle,72)
        assert rsvg.rsvg_handle_render_document(handle,context,ct.byref(Rectangle(0,0,width,height)),None)
        cairo.cairo_surface_finish(surface)
        assert cairo.cairo_surface_status(surface)==0
    finally:
        cairo.cairo_destroy(context)
        cairo.cairo_surface_destroy(surface)
        gobject.g_object_unref(handle)
