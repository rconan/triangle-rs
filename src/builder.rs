use std::{ffi::CString, fmt, slice};

use crate::{
    ffi::{triangulate, triangulateio, trifree},
    Delaunay,
};

/// Delaunay triangulation builder
#[derive(Debug)]
#[allow(dead_code)]
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
    pub fn add_nodes(mut self, nodes: &[f64]) -> Self {
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
    pub fn set_switches(mut self, switches: &str) -> Self {
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
        //        println!("Delaunay triangulation with  switches: {}", self.switches);
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
            let buffer = slice::from_raw_parts(delaunay.pointlist, n).to_vec();
            trifree(delaunay.pointlist as _);
            buffer
        };
        let point_markers: Vec<i32> = unsafe {
            let n = delaunay.numberofpoints as usize;
            let buffer = slice::from_raw_parts(delaunay.pointmarkerlist, n).to_vec();
            trifree(delaunay.pointmarkerlist as _);
            buffer
        };
        let triangles: Vec<usize> = unsafe {
            let n = delaunay.numberoftriangles as usize * 3;
            let buffer = slice::from_raw_parts(delaunay.trianglelist, n).to_vec();
            trifree(delaunay.trianglelist as _);
            buffer
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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn area() {
        let tri = Builder::new()
            .add_nodes(&[0., 0.])
            .add_polygon(&[1., 0., 0., 1., -1., 0., 0., -1.])
            .build();
        assert_eq!(tri.area(), 2.)
    }
}
