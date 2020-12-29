use plotters::prelude::*;
use std::ffi::CString;
use std::path::Path;

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
    pub triangles: Vec<usize>,
}

pub trait TriDraw {
    fn mesh<T: AsRef<Path>>(&self, path: T, lim: f64);
}
impl TriDraw for Delaunay {
    fn mesh<T: AsRef<Path>>(&self, path: T, lim: f64) {
        let p_x: Vec<_> = self.points.chunks(2).map(|x| x[0]).collect();
        let p_y: Vec<_> = self.points.chunks(2).map(|x| x[1]).collect();
        let plot = SVGBackend::new(&path, (768, 768)).into_drawing_area();
        plot.fill(&WHITE).unwrap();
        let mut chart = ChartBuilder::on(&plot)
            .set_label_area_size(LabelAreaPosition::Left, 40)
            .set_label_area_size(LabelAreaPosition::Bottom, 40)
            .build_cartesian_2d(-lim..lim, -lim..lim)
            .unwrap();
        chart.configure_mesh().draw().unwrap();
        let vertices: Vec<Vec<(f64, f64)>> = self
            .triangles
            .chunks(3)
            .map(|t| {
                t.iter()
                    .map(|&i| (p_x[i as usize], p_y[i as usize]))
                    .collect::<Vec<(f64, f64)>>()
            })
            .collect();
        /*
        vertices.iter().for_each(|v| {
            chart
                .draw_series(std::iter::once(Polygon::new(v.clone(), &RED.mix(0.2))))
                .unwrap();
        });
        */
        vertices.iter().for_each(|v| {
            chart
                .draw_series(LineSeries::new(
                    v.iter().cycle().take(4).map(|(x, y)| (*x, *y)),
                    &BLACK,
                ))
                .unwrap();
        });
        /*
        chart
            .draw_series(
                p_x.iter()
                    .cloned()
                    .zip(p_y.iter().cloned())
                    .map(|p| Circle::new(p, 3, RED.filled())),
            )
            .unwrap();
        */
    }
}

impl Delaunay {
    pub fn new() -> Self {
        Self {
            points: vec![],
            triangles: vec![],
        }
    }
    pub fn n_triangles(&self) -> usize {
        self.triangles.len() / 3
    }
    pub fn is_point_inside(&self, point: &[f64], triangle_id: usize) -> bool {
        let triangle = self.triangles.chunks(3).nth(triangle_id).unwrap();
        let points: Vec<&[f64]> = self.points.chunks(2).collect();
        for i in 0..3 {
            let j = (i + 1) % 3;
            let vi = triangle[i];
            let vj = triangle[j];
            let d = (points[vj][0] - points[vi][0]) * (point[1] - points[vi][1])
                - (points[vj][1] - points[vi][1]) * (point[0] - points[vi][0]);
            if d < 0. && d.abs() > 1e-9 {
                return false;
            }
        }
        true
    }
    pub fn which_contains_point(&self, point: &[f64]) -> Option<usize> {
        for k in 0..self.n_triangles() {
            if self.is_point_inside(point, k) {
                return Some(k);
            }
        }
        None
    }
    pub fn point_into_barycentric(&self, point: &[f64], triangle_ids: &[usize]) -> [f64; 3] {
        let points: Vec<&[f64]> = self.points.chunks(2).collect();
        let v: Vec<&[f64]> = triangle_ids.iter().map(|&i| points[i]).collect();
        let area =
            (v[1][1] - v[2][1]) * (v[0][0] - v[2][0]) + (v[2][0] - v[1][0]) * (v[0][1] - v[2][1]);
        let w0 = ((v[1][1] - v[2][1]) * (point[0] - v[2][0])
            + (v[2][0] - v[1][0]) * (point[1] - v[2][1]))
            / area;
        let w1 = ((v[2][1] - v[0][1]) * (point[0] - v[2][0])
            + (v[0][0] - v[2][0]) * (point[1] - v[2][1]))
            / area;
        [w0, w1, 1. - w0 - w1]
    }
    pub fn barycentric_interpolation(&self, point: &[f64], val_at_points: &[f64]) -> f64 {
        match self.which_contains_point(point) {
            Some(tid) => {
                let triangle_ids = self.triangles.chunks(3).nth(tid).unwrap();
                let values: Vec<f64> = triangle_ids.iter().map(|&i| val_at_points[i]).collect();
                self.point_into_barycentric(point, triangle_ids)
                    .iter()
                    .zip(values.iter())
                    .fold(0., |a, (w, v)| a + w * v)
            }
            None => std::f64::NAN,
        }
    }
    /*
    pub fn linear_interpolation(&self, point: &[f64], val_at_points: &[f64]) -> f64 {
        match self.contain_point(point) {
            Some(tid) => {
                println!("Triangle #{}", tid);
                let ipts = self.triangles[tid].clone();
                let mut wgts = vec![0f64; 3];
                for i in 0..3 {
                    //let j = if i + 1 == 3 { 0 } else { i + 1 };
                    //let k = if i + 2 == 3 { 0 } else { i + 2 };
                    let j = (i + 1) % 3;
                    let k = (i + 2) % 3;
                    wgts[k] = (self.points[ipts[j]][0] - self.points[ipts[i]][0])
                        * (point[1] - self.points[ipts[i]][1])
                        - (self.points[ipts[j]][1] - self.points[ipts[i]][1])
                            * (point[0] - self.points[ipts[i]][0]);
                }
                println!("weights: {:?}", wgts);
                let sum: f64 = wgts.iter().sum();
                let values: Vec<f64> = ipts.iter().map(|k| val_at_points[*k]).collect();
                println!("values: {:?}", values);
                wgts.iter()
                    .zip(values.iter())
                    .fold(0f64, |a, (w, v)| a + w * v / sum)
            }
            None => std::f64::NAN,
        }
    }
    */
}

pub enum TriangulateIO {
    Points(Vec<f64>),
}

pub struct Builder {
    triangulate_io: Vec<TriangulateIO>,
    switches: String,
}
impl Builder {
    pub fn new() -> Self {
        Self {
            triangulate_io: vec![],
            switches: "z".to_owned(),
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
            ..self
        }
    }
    pub fn set_switches(self, switches: &str) -> Self {
        Self {
            switches: format!("z{}", switches),
            ..self
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
                    tri_io.numberofpoints = p.len() as i32 / 2;
                    tri_io.pointlist = p.as_mut_ptr();
                }
            }
        }
        let mut delaunay: triangulateio = unsafe { std::mem::zeroed() };
        let switches = CString::new(self.switches).unwrap();
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
            let n = delaunay.numberofpoints as usize * 2;
            Vec::from_raw_parts(delaunay.pointlist, n, n)
        };
        let triangles: Vec<usize> = unsafe {
            let n = delaunay.numberoftriangles as usize * 3;
            Vec::from_raw_parts(delaunay.trianglelist, n, n)
        }
        .iter()
        .map(|x| *x as usize)
        .collect();
        Delaunay { points, triangles }
    }
}
