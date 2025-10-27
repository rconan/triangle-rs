#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct triangulateio {
    pub pointlist: *mut f64,
    pub pointattributelist: *mut f64,
    pub pointmarkerlist: *mut ::std::os::raw::c_int,
    pub numberofpoints: ::std::os::raw::c_int,
    pub numberofpointattributes: ::std::os::raw::c_int,
    pub trianglelist: *mut ::std::os::raw::c_int,
    pub triangleattributelist: *mut f64,
    pub trianglearealist: *mut f64,
    pub neighborlist: *mut ::std::os::raw::c_int,
    pub numberoftriangles: ::std::os::raw::c_int,
    pub numberofcorners: ::std::os::raw::c_int,
    pub numberoftriangleattributes: ::std::os::raw::c_int,
    pub segmentlist: *mut ::std::os::raw::c_int,
    pub segmentmarkerlist: *mut ::std::os::raw::c_int,
    pub numberofsegments: ::std::os::raw::c_int,
    pub holelist: *mut f64,
    pub numberofholes: ::std::os::raw::c_int,
    pub regionlist: *mut f64,
    pub numberofregions: ::std::os::raw::c_int,
    pub edgelist: *mut ::std::os::raw::c_int,
    pub edgemarkerlist: *mut ::std::os::raw::c_int,
    pub normlist: *mut f64,
    pub numberofedges: ::std::os::raw::c_int,
}
impl Default for triangulateio {
    fn default() -> Self {
        triangulateio {
            pointlist: std::ptr::null_mut::<f64>(),
            pointattributelist: std::ptr::null_mut::<f64>(),
            pointmarkerlist: std::ptr::null_mut::<i32>(),
            numberofpoints: 0i32,
            numberofpointattributes: 0i32,
            trianglelist: std::ptr::null_mut::<i32>(),
            triangleattributelist: std::ptr::null_mut::<f64>(),
            trianglearealist: std::ptr::null_mut::<f64>(),
            neighborlist: std::ptr::null_mut::<i32>(),
            numberoftriangles: 0i32,
            numberofcorners: 0i32,
            numberoftriangleattributes: 0i32,
            segmentlist: std::ptr::null_mut::<i32>(),
            segmentmarkerlist: std::ptr::null_mut::<i32>(),
            numberofsegments: 0i32,
            holelist: std::ptr::null_mut::<f64>(),
            numberofholes: 0i32,
            regionlist: std::ptr::null_mut::<f64>(),
            numberofregions: 0i32,
            edgelist: std::ptr::null_mut::<i32>(),
            edgemarkerlist: std::ptr::null_mut::<i32>(),
            normlist: std::ptr::null_mut::<f64>(),
            numberofedges: 0i32,
        }
    }
}

extern "C" {
    pub fn triangulate(
        arg1: *mut ::std::os::raw::c_char,
        arg2: *mut triangulateio,
        arg3: *mut triangulateio,
        arg4: *mut triangulateio,
    );
}
extern "C" {
    pub fn trifree(memptr: *mut ::std::os::raw::c_int);
}
