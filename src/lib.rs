use std::ffi::CString;

//include!("bindings.rs");

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

#[derive(Debug)]
pub struct Delaunay {
    pub points: Vec<f64>,
    pub triangles: Vec<Vec<i32>>,
}
impl Delaunay {
    pub fn new() -> Self {
        Self {
            points: vec![],
            triangles: vec![],
        }
    }
}

pub enum TriangulateIO {
    Points(Vec<f64>),
}

pub struct Builder {
    triangulate_io: Vec<TriangulateIO>,
}

impl Builder {
    pub fn new() -> Self {
        Self {
            triangulate_io: vec![],
        }
    }
    pub fn set_points(self, x: Vec<f64>, y: Vec<f64>) -> Self {
        assert!(x.len() == y.len(), "x and y are not the same length.");
        let mut data = self.triangulate_io;
        let xy = x
            .into_iter()
            .zip(y.into_iter())
            .flat_map(|(x, y)| vec![x, y])
            .collect();
        data.push(TriangulateIO::Points(xy));
        Self {
            triangulate_io: data,
        }
    }
    pub fn build(self) -> Delaunay {
        use TriangulateIO::*;
        let mut tri_io: triangulateio = unsafe { std::mem::zeroed() };
        tri_io.numberofpoints = 0_i32;
        tri_io.numberofpointattributes = 0_i32;
        tri_io.numberoftriangles = 0_i32;
        tri_io.numberofcorners = 0_i32;
        tri_io.numberoftriangleattributes = 0_i32;
        tri_io.numberofsegments = 0_i32;
        tri_io.numberofholes = 0_i32;
        tri_io.numberofregions = 0_i32;
        tri_io.numberofedges = 0_i32;
        for t in self.triangulate_io {
            match t {
                Points(mut p) => {
                    println!("p: {:?}", p);
                    tri_io.numberofpoints = p.len() as i32 / 2;
                    tri_io.pointlist = p.as_mut_ptr();
                }
            }
        }
        let mut delaunay: triangulateio = unsafe { std::mem::zeroed() };
        let mut dtri = Delaunay::new();
        //delaunay.pointlist = dtri.points.as_mut_ptr() as *mut _;
        let switches = CString::new("z").unwrap();
        unsafe {
            let mut empty_tri: triangulateio = std::mem::zeroed();
            triangulate(
                switches.into_raw(),
                &mut tri_io,
                &mut delaunay,
                &mut empty_tri,
            )
        };
        let points: Vec<f64> = unsafe {
            let n = delaunay.numberofpoints as usize;
            Vec::from_raw_parts(delaunay.pointlist, n, n)
        };
        let triangles: Vec<Vec<i32>> = unsafe {
            let n = delaunay.numberoftriangles as usize * 3;
            Vec::from_raw_parts(delaunay.trianglelist, n, n)
        }
        .chunks(3)
        .map(|x| x.to_vec())
        .collect();
        Delaunay { points, triangles }
    }
}
