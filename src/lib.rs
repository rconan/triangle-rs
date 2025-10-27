//! A Rust library wrapped around the 2D mesh generator and Delaunay triangulator [Triangle](https://www.cs.cmu.edu/~quake/triangle.html)

use rayon::prelude::*;
use serde::Serialize;
use serde_pickle as pkl;
use std::{error::Error, fmt, fs::File, path::Path};

mod builder;
mod ffi;
pub use builder::Builder;

/// Delaunay triangulation
#[derive(Default, Debug, Serialize)]
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
    /// Returns the `Delaunay` builder
    pub fn builder() -> Builder {
        Builder::new()
    }
    /// Returns the number of vertices
    pub fn n_vertices(&self) -> usize {
        self.points.len() / 2
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
    /// Returns an interator over the triangles, each item is a vector of the 3 (x,y) vertices coordinates
    pub fn triangle_vertex_iter(&self) -> impl Iterator<Item = Vec<(f64, f64)>> + '_ {
        let x = self.x();
        let y = self.y();
        self.triangle_iter()
            .map(move |t| t.iter().map(|&i| (x[i], y[i])).collect::<Vec<(f64, f64)>>())
    }

    /// Returns the triangle areas
    pub fn triangle_areas(&self) -> Vec<f64> {
        let vertices: Vec<Vec<f64>> = self.vertex_iter().map(|x| x.to_vec()).collect();
        self.triangle_iter()
            .map(|t| {
                let (a, b, c) = (&vertices[t[0]], &vertices[t[1]], &vertices[t[2]]);
                0.5 * ((a[0] - c[0]) * (b[1] - a[1]) - (a[0] - b[0]) * (c[1] - a[1])).abs()
            })
            .collect()
    }
    /// Returns the area covered by the mesh as the sum of the triangle area
    pub fn area(&self) -> f64 {
        let vertices: Vec<Vec<f64>> = self.vertex_iter().map(|x| x.to_vec()).collect();
        self.triangle_iter().fold(0., |s, t| {
            let (a, b, c) = (&vertices[t[0]], &vertices[t[1]], &vertices[t[2]]);
            s + 0.5 * ((a[0] - c[0]) * (b[1] - a[1]) - (a[0] - b[0]) * (c[1] - a[1])).abs()
        })
    }
    /// Returns the area covered by the mesh as the sum of the Delaunay triangles area
    pub fn mesh_area(&self) -> f64 {
        let vertices: Vec<Vec<f64>> = self.vertex_iter().map(|x| x.to_vec()).collect();
        self.triangle_iter().fold(0., |s, t| {
            let (a, b, c) = (&vertices[t[0]], &vertices[t[1]], &vertices[t[2]]);
            s + 0.5 * ((a[0] - c[0]) * (b[1] - a[1]) - (a[0] - b[0]) * (c[1] - a[1])).abs()
        })
    }
    pub fn average(&self, vertices: &[[f64; 3]], data: &[f64]) -> f64 {
        self.triangle_iter().fold(0., |s, t| {
            let (a, b, c) = (&vertices[t[0]], &vertices[t[1]], &vertices[t[2]]);
            let ta = 0.5 * ((a[0] - c[0]) * (b[1] - a[1]) - (a[0] - b[0]) * (c[1] - a[1])).abs();
            let sa = t.iter().fold(0., |m, i| m + data[*i]) / 3 as f64;
            s + ta * sa
        })
    }
    pub fn average_with<F: Fn(f64) -> f64>(
        &self,
        vertices: &[[f64; 3]],
        data: &[f64],
        f: F,
    ) -> f64 {
        self.triangle_iter().fold(0., |s, t| {
            let (a, b, c) = (&vertices[t[0]], &vertices[t[1]], &vertices[t[2]]);
            let ta = 0.5 * ((a[0] - c[0]) * (b[1] - a[1]) - (a[0] - b[0]) * (c[1] - a[1])).abs();
            let sa = t.iter().fold(0., |m, i| m + f(data[*i])) / 3 as f64;
            s + ta * sa
        })
    }
    /// Returns the dot product of two vector `a` and `b` sampled on the Delaunay mesh
    pub fn dot(&self, v_a: &[f64], v_b: &[f64]) -> f64 {
        assert_eq!(v_a.len(), self.n_vertices());
        assert_eq!(v_b.len(), self.n_vertices());
        let vertices: Vec<Vec<f64>> = self.vertex_iter().map(|x| x.to_vec()).collect();
        self.triangle_iter().fold(0., |s, t| {
            let (a, b, c) = (&vertices[t[0]], &vertices[t[1]], &vertices[t[2]]);
            let ta = 0.5 * ((a[0] - c[0]) * (b[1] - a[1]) - (a[0] - b[0]) * (c[1] - a[1])).abs();
            let sa = t.iter().fold(0., |m, i| m + v_a[*i]) / 3 as f64;
            let sb = t.iter().fold(0., |m, i| m + v_b[*i]) / 3 as f64;
            s + ta * sa * sb
        }) / self.area()
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
    pub fn dump<T: AsRef<Path>>(&self, filename: T) -> Result<(), Box<dyn Error>> {
        pkl::to_writer(&mut File::create(filename)?, self, true)?;
        Ok(())
    }
}

impl fmt::Display for Delaunay {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let areas = self.triangle_areas();
        let areas_min = areas.iter().cloned().fold(f64::INFINITY, f64::min);
        let areas_max = areas.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let areas_sum = areas.iter().sum::<f64>();
        let areas_mean = areas_sum / areas.len() as f64;
        write!(
            f,
            r#"Delaunay triangulation:
 - vertices: {}
 - triangles: {}
 - triangles area: 
   - min : {:.6}
   - max : {:.6}
   - mean: {:.6}
   - sum : {:.6}"#,
            self.n_vertices(),
            self.n_triangles(),
            areas_min,
            areas_max,
            areas_mean,
            areas_sum
        )
    }
}

pub enum TriangulateIO {
    Points(Vec<f64>),
}

/*
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
    fn heatmap(
        &self,
        _x: &[f64],
        _y: &[f64],
        _z: &[f64],
        _range: Range<f64>,
        _config: Option<Config>,
    ) -> Result<(), Box<dyn Error>> {
        Ok(())
    }
}
 */

