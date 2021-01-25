//! A Rust library wrapped around the 2D mesh generator and Delaunay triangulator [Triangle](https://www.cs.cmu.edu/~quake/triangle.html)

use complot::TriPlot;
use plotters::coord::types::RangedCoordf64;
use plotters::prelude::{Cartesian2d, ChartContext, DrawingBackend, LineSeries, RGBColor};
use rayon::prelude::*;
use serde::Serialize;
use std::ffi::CString;
use std::fmt;

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

/// Delaunay triangulation
#[derive(Debug, Serialize)]
pub struct Delaunay {
    /// Triangulation vertices as [x0,y0,x1,y1,...]
    pub points: Vec<f64>,
    /// Indices in `points.chunks(2)`, the first 3 indices correspond to the vertices of the 1st triangle, the next 3 to the 2nd triangle, etc...
    pub point_markers: Vec<i32>,
    pub triangles: Vec<usize>,
    /// List of triangles neighbors: indices in triangles.chunks(3) (3 integers per triangle)
    pub neighbors: Option<Vec<i32>>,
    /// Edges endpoints: indices in points.chunks(2) (2 integers per edge)
    pub edges: Option<Vec<usize>>,
}

impl Delaunay {
    /// Creates a new empty Delaunay triangulation
    pub fn new() -> Self {
        Self {
            points: vec![],
            point_markers: vec![],
            triangles: vec![],
            neighbors: None,
            edges: None,
        }
    }
    /// Returns the number of Delaunay triangles
    pub fn n_triangles(&self) -> usize {
        self.triangles.len() / 3
    }
    /// Returns an iterator over the vertices, each item is a vertex (x,y) coordinates
    pub fn vertex_iter(&self) -> std::slice::Chunks<'_, f64> {
        self.points.chunks(2)
    }
    /// Returns an iterator over mutable vertices, each item is a vertex (x,y) coordinates
    pub fn vertex_iter_mut(&mut self) -> std::slice::ChunksMut<'_, f64> {
        self.points.chunks_mut(2)
    }
    /// Returns a parallel iterator over the vertices, each item is a vertex (x,y) coordinates
    pub fn vertex_par_iter(&self) -> rayon::slice::Chunks<'_, f64> {
        self.points.par_chunks(2)
    }
    /// Returns an iterator over the triangles, each item is the indices of the vertices in `vertex_iter`
    pub fn triangle_iter(&self) -> std::slice::Chunks<'_, usize> {
        self.triangles.chunks(3)
    }
    /// Returns an iterator over mutable triangles, each item is the indices of the vertices in `vertex_iter`
    pub fn triangle_iter_mut(&mut self) -> std::slice::ChunksMut<'_, usize> {
        self.triangles.chunks_mut(3)
    }
    /// Returns a parallel iterator over the triangles, each item is the indices of the vertices in `vertex_iter`
    pub fn triangle_par_iter(&self) -> rayon::slice::Chunks<'_, usize> {
        self.triangles.par_chunks(3)
    }
    /// Gets node x coordinates
    pub fn x(&self) -> Vec<f64> {
        self.vertex_iter().map(|xy| xy[0]).collect()
    }
    /// Gets node y coordinates
    pub fn y(&self) -> Vec<f64> {
        self.vertex_iter().map(|xy| xy[1]).collect()
    }
    // Returns the area covered by the mesh as the sum of the triangle area
    pub fn area(&self) -> f64 {
        let vertices: Vec<Vec<f64>> = self.vertex_iter().map(|x| x.to_vec()).collect();
        self.triangle_iter().fold(0., |s, t| {
            let (a, b, c) = (&vertices[t[0]], &vertices[t[1]], &vertices[t[2]]);
            s + 0.5 * ((a[0] - c[0]) * (b[1] - a[1]) - (a[0] - b[0]) * (c[1] - a[1])).abs()
        })
    }
    /// Returns true if a point `[x,y]` is inside the triangle given by its index (`triangle_id`) in `triangles_iter`, otherwise returns false
    pub fn is_point_inside(&self, point: &[f64], triangle_id: usize) -> bool {
        let triangle = self.triangle_iter().nth(triangle_id).unwrap();
        let points: Vec<&[f64]> = self.vertex_iter().collect();
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
    /// Finds the index of the triangle in `triangles_iter` that contains the given point `[x,y]`
    pub fn which_contains_point(&self, point: &[f64]) -> Option<usize> {
        for k in 0..self.n_triangles() {
            if self.is_point_inside(point, k) {
                return Some(k);
            }
        }
        None
    }
    /// Returns the barycentric coordinates of a point `[x,y]` with respect to the triangle that contains it
    ///
    /// The triangle that contains the point is specified with the indices `triangle_ids` in `vertex_iter` of the triangle vertices
    pub fn point_into_barycentric(&self, point: &[f64], triangle_ids: &[usize]) -> [f64; 3] {
        let points: Vec<&[f64]> = self.vertex_iter().collect();
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
    /// Linearly interpolates at a given point [x,y], values `val_at_vertices` at the Delaunay mesh vertices
    ///
    /// The linear interpolation is based on the barycentric coordinates of the point
    pub fn barycentric_interpolation(&self, point: &[f64], val_at_vertices: &[f64]) -> f64 {
        match self.which_contains_point(point) {
            Some(tid) => {
                let triangle_ids = self.triangle_iter().nth(tid).unwrap();
                let values: Vec<f64> = triangle_ids.iter().map(|&i| val_at_vertices[i]).collect();
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

/// Delaunay triangulation builder
#[derive(Debug)]
pub struct Builder {
    //triangulate_io: Vec<TriangulateIO>,
    points: Vec<f64>,
    segments: Option<Vec<i32>>,
    n_segments: i32,
    holes: Option<Vec<f64>>,
    n_holes: i32,
    switches: String,
    boundary_marker: i32,
    point_markers: Option<Vec<i32>>,
    segment_markers: Option<Vec<i32>>,
    tri_io: triangulateio,
}
impl Builder {
    /// Creates a new Delaunay triangulation builder
    pub fn new() -> Self {
        Self {
            //triangulate_io: vec![],
            switches: "z".to_owned(),
            points: vec![],
            segments: None,
            n_segments: 032,
            holes: None,
            n_holes: 0i32,
            boundary_marker: 1i32,
            point_markers: None,
            segment_markers: None,
            tri_io: triangulateio::default(),
        }
    }
    /// Sets the Delaunay mesh `x` and `y` vertices coordinates
    pub fn add_nodes(&mut self, nodes: &[f64]) -> &mut Self {
        self.points.extend(nodes);
        self
    }
    pub fn set_segments(self, x: Vec<i32>, y: Vec<i32>) -> Self {
        assert!(x.len() == y.len(), "x and y are not the same length.");
        //let mut data = self.triangulate_io;
        let n = x.len() as i32;
        let xy = x
            .into_iter()
            .zip(y.into_iter())
            .flat_map(|(x, y)| vec![x, y])
            .collect();
        //data.push(TriangulateIO::Points(xy));
        Self {
            segments: Some(xy),
            n_segments: n,
            ..self
        }
    }
    pub fn add_holes(&mut self, x: f64, y: f64) -> &mut Self {
        match self.holes {
            Some(ref mut h) => {
                h.extend(vec![x, y]);
            }
            None => {
                self.holes = Some(vec![x, y]);
            }
        }
        self
    }
    /// Sets the Delaunay mesh vertices as [x0,y0,x1,y1,...]
    pub fn set_tri_points(self, points: Vec<f64>) -> Self {
        /*let mut data = self.triangulate_io;
        data.push(TriangulateIO::Points(points));*/
        Self { points, ..self }
    }
    /// Adds a closed polygon given its vertices [x1,y1,x2,y2,...]
    pub fn add_polygon(&mut self, vertices: &[f64]) -> &mut Self {
        //let boundary_marker = self.boundary_marker + 1;
        let a = (self.points.len() / 2) as i32;
        let n_segments = (vertices.len() / 2) as i32;
        /*
        let point_markers = match self.point_markers.clone() {
            Some(mut p_m) => {
                p_m.extend(vec![boundary_marker; n_segments as usize]);
                p_m
            }
            None => vec![boundary_marker; n_segments as usize],
        };
        let segment_markers = match self.segment_markers.clone() {
            Some(mut s_m) => {
                s_m.extend(vec![boundary_marker; n_segments as usize]);
                s_m
            }
            None => vec![boundary_marker; n_segments as usize],
        };
        println!("point markers: {:?}", point_markers);*/
        let segments_vertices = (0..n_segments).flat_map(|k| vec![a + k, a + (k + 1) % n_segments]);
        match self.segments {
            Some(ref mut s) => {
                s.extend(segments_vertices);
            }
            None => {
                self.segments = Some(segments_vertices.collect::<Vec<i32>>());
            }
        };
        self.points.extend(vertices);
        self
    }
    /// Sets triangulation [switches](https://www.cs.cmu.edu/~quake/triangle.switch.html)
    pub fn set_switches(&mut self, switches: &str) -> &mut Self {
        self.switches = format!("z{}", switches);
        self
    }
    /// Compute the Delaunay mesh and returns a `Delaunay` structure
    pub fn build(&mut self) -> Delaunay {
        self.tri_io.numberofpoints = (self.points.len() / 2) as i32;
        self.tri_io.pointlist = self.points.as_mut_ptr();
        if let Some(ref mut s) = self.segments {
            self.tri_io.numberofsegments = (s.len() / 2) as i32;
            self.tri_io.segmentlist = s.as_mut_ptr();
        }
        if let Some(ref mut h) = self.holes {
            self.tri_io.numberofholes = (h.len() / 2) as i32;
            self.tri_io.holelist = h.as_mut_ptr();
        }
        //use TriangulateIO::*;
        let mut delaunay: triangulateio = unsafe { std::mem::zeroed() };
        let switches = CString::new(self.switches.as_str()).unwrap();
        unsafe {
            let mut empty_tri: triangulateio = std::mem::zeroed();
            triangulate(
                switches.into_raw(),
                &mut self.tri_io,
                &mut delaunay,
                &mut empty_tri,
            )
        };
        let points: Vec<f64> = unsafe {
            let n = delaunay.numberofpoints as usize * 2;
            Vec::from_raw_parts(delaunay.pointlist, n, n)
        };
        let point_markers: Vec<i32> = unsafe {
            let n = delaunay.numberofpoints as usize;
            Vec::from_raw_parts(delaunay.pointmarkerlist, n, n)
        };
        let triangles: Vec<usize> = unsafe {
            let n = delaunay.numberoftriangles as usize * 3;
            Vec::from_raw_parts(delaunay.trianglelist, n, n)
        }
        .iter()
        .map(|x| *x as usize)
        .collect();
        let neighbors: Option<Vec<i32>> = if self.switches.contains("n") {
            let n = delaunay.numberoftriangles as usize * 3;
            Some(unsafe { Vec::from_raw_parts(delaunay.neighborlist, n, n) })
        } else {
            None
        };
        let edges: Option<Vec<usize>> = if self.switches.contains("e") {
            let n = delaunay.numberofedges as usize * 2;
            Some(
                unsafe { Vec::from_raw_parts(delaunay.edgelist, n, n) }
                    .iter()
                    .map(|x| *x as usize)
                    .collect(),
            )
        } else {
            None
        };
        Delaunay {
            points,
            point_markers,
            triangles,
            neighbors,
            edges,
        }
    }
}
impl From<Vec<f64>> for Builder {
    fn from(points: Vec<f64>) -> Self {
        Self {
            //triangulate_io: vec![TriangulateIO::Points(points)],
            points,
            switches: "z".to_owned(),
            segments: None,
            n_segments: 032,
            holes: None,
            n_holes: 0i32,
            boundary_marker: 1i32,
            point_markers: None,
            segment_markers: None,
            tri_io: triangulateio::default(),
        }
    }
}
impl From<&[f64]> for Builder {
    fn from(points: &[f64]) -> Self {
        Self {
            //triangulate_io: vec![TriangulateIO::Points(points.to_owned())],
            points: points.to_owned(),
            switches: "z".to_owned(),
            segments: None,
            n_segments: 032,
            holes: None,
            n_holes: 0i32,
            boundary_marker: 1i32,
            point_markers: None,
            segment_markers: None,
            tri_io: triangulateio::default(),
        }
    }
}

impl fmt::Display for Builder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let nodes_1st_line = format!("{} 2 0 0", self.points.len() / 2);
        let nodes = self
            .points
            .chunks(2)
            .enumerate()
            .map(|(k, xy)| format!("{}  {} {}", k, xy[0], xy[1]))
            .collect::<Vec<String>>()
            .join("\n");
        let segs = match &self.segments {
            Some(s) => {
                let segs_1st_line = format!("{} 0", s.len() / 2);
                let segs = s
                    .chunks(2)
                    .enumerate()
                    .map(|(k, xy)| format!("{}  {} {}", k, xy[0], xy[1]))
                    .collect::<Vec<String>>()
                    .join("\n");
                let holes_1st_line = format!("{}", self.n_holes);
                let holes = match &self.holes {
                    Some(h) => h
                        .chunks(2)
                        .enumerate()
                        .map(|(k, xy)| format!("{}  {} {}", k, xy[0], xy[1]))
                        .collect::<Vec<String>>()
                        .join("\n"),
                    None => "".to_owned(),
                };
                [segs_1st_line, segs, holes_1st_line, holes].join("\n")
            }
            None => "".to_owned(),
        };
        write!(f, "{}", [nodes_1st_line, nodes, segs].join("\n"))
    }
}

impl TriPlot for Delaunay {
    fn mesh<'a, D: DrawingBackend>(
        &self,
        x: &[f64],
        y: &[f64],
        color: [u8; 3],
        chart: &mut ChartContext<'a, D, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    ) -> &Self {
        let color = RGBColor(color[0], color[1], color[2]);
        self.triangle_iter()
            .map(|t| t.iter().map(|&i| (x[i], y[i])).collect::<Vec<(f64, f64)>>())
            .for_each(|v| {
                chart
                    .draw_series(LineSeries::new(
                        v.iter().cycle().take(4).map(|(x, y)| (*x, *y)),
                        &color,
                    ))
                    .unwrap();
            });
        self
    }
    fn map<'a, D: DrawingBackend>(
        &self,
        _x: &[f64],
        _y: &[f64],
        _z: &[f64],
        _chart: &mut ChartContext<'a, D, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    ) -> &Self {
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn area() {
        let tri = Builder::new().add_nodes(&[0.,0.]).add_polygon(&[1.,0.,0.,1.,-1.,0.,0.,-1.]).build();
        println!("area: {}",tri.area());
    }
}
