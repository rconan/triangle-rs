use triangle_rs::Builder;

#[derive(Clone)]
struct TrussVertices {
    x: Vec<f64>,
    y: Vec<f64>,
}
impl TrussVertices {
    fn zip(&self) -> std::iter::Zip<std::slice::Iter<'_, f64>, std::slice::Iter<'_, f64>> {
        self.x.iter().zip(self.y.iter())
    }
    fn into_zip(self) -> std::iter::Zip<std::vec::IntoIter<f64>, std::vec::IntoIter<f64>> {
        self.x.into_iter().zip(self.y.into_iter())
    }
    fn zip_mut(
        &mut self,
    ) -> std::iter::Zip<std::slice::IterMut<'_, f64>, std::slice::IterMut<'_, f64>> {
        self.x.iter_mut().zip(self.y.iter_mut())
    }
}

fn main() {
    let d = 0.25;
    let outer_radius = 8.365 / 2.;
    let inner_radius = 1.675;

    let mut truss_vertices = TrussVertices {
        x: vec![
            -3.011774, -2.446105, -3.011774, -2.799304, -2.33903, -1.566412, -1.640648, -1.9445,
            -1.640648, -1.566412, -2.347462, -1.597649, -1.725044, -2.392888, -2.799304,
        ],
        y: vec![
            -2.902158, 0., 2.902158, 3.107604, 0.07244, 0.518512, 0.175429, 0., -0.175429,
            -0.518512, -0.067572, -3.865336, -3.810188, -0.427592, -3.107604,
        ],
    };
    truss_vertices.zip_mut().for_each(|(x, y)| {
        let r = x.hypot(*y);
        if r < inner_radius {
            let (s, c) = y.atan2(*x).sin_cos();
            *x = r * c;
            *y = r * s;
        }
    });

    let trusses: Vec<TrussVertices> = (0..3)
        .map(|k| {
            let (s, c) = ((k + 1) as f64 * (-2. * std::f64::consts::PI / 3.)).sin_cos();
            let mut truss = TrussVertices {
                x: vec![],
                y: vec![],
            };
            truss_vertices.zip().for_each(|(x, y)| {
                let xp = x * c + y * s;
                let yp = -x * s + y * c;
                truss.x.push(xp);
                truss.y.push(yp);
            });
            truss
        })
        .collect();

    let xy_arc = (0..3).flat_map(|k| {
        print!("Truss pair: [{},{}], ", k, (k + 1) % 3);
        let t_a = &trusses[k];
        let t_b = &trusses[(k + 1) % 3];
        let o1 = t_a.y[9].atan2(t_a.x[9]);
        let o2 = t_b.y[5].atan2(t_b.x[5]);
        let n = (inner_radius * (o2 - o1) / d).ceil() as u32;
        print!(
            "arc: [{:+7.2}, {:+7.2}, {}], ",
            o1.to_degrees(),
            o2.to_degrees(),
            n
        );
        let mut xy_arc: Vec<_> = (0..n)
            .map(|k| {
                let (s, c) = (o1 + (o2 - o1) * (k as f64 / ((n - 1) as f64))).sin_cos();
                (inner_radius * c, inner_radius * s)
            })
            .collect();
        let o1 = t_b.y[6].atan2(t_b.x[6]);
        let o2 = {
            let o = t_b.y[8].atan2(t_b.x[8]);
            if o < o1 {
                2. * std::f64::consts::PI + o
            } else {
                o
            }
        };
        let n = (inner_radius * (o2 - o1) / d).ceil() as i32;
        println!(
            "arc: [{:+7.2}, {:+7.2}, {}]",
            o1.to_degrees(),
            o2.to_degrees(),
            n
        );
        xy_arc.extend((0..n).map(|k| {
            let (s, c) = (o1 + (o2 - o1) * (k as f64 / ((n - 1) as f64))).sin_cos();
            (inner_radius * c, inner_radius * s)
        }));
        xy_arc
    });

    // NODES
    let n = (2. * std::f64::consts::PI * outer_radius / d).ceil() as usize;
    let mut nodes: Vec<_> = (0..n)
        .map(|k| {
            let (s, c) = (2. * std::f64::consts::PI * k as f64 / n as f64).sin_cos();
            (outer_radius * c, outer_radius * s)
        })
        .collect();
    nodes.extend(trusses.iter().flat_map(|t| t.zip().map(|(x, y)| (*x, *y))));
    nodes.extend(xy_arc.clone());
    // EDGES
    let mut edges: Vec<_> = (0..n).map(|k| (k as i32, ((k + 1) % n) as i32)).collect();
    edges.extend(
        trusses
            .iter()
            .flat_map(|t| t.clone().into_zip())
            .chain(xy_arc.collect::<Vec<(f64, f64)>>())
            .enumerate()
            .map(|(k, _)| ((n + k) as i32, (n + k + 1) as i32)),
    );
    // HOLES
    let holes = [[-2.2, 0.], [1.1, -1.9], [1.1, 1.9], [0., 0.]];

    println!("nodes: {:#?}",nodes);
    let (xn, yn): (Vec<f64>, Vec<f64>) = nodes.into_iter().unzip();
    let (xe, ye): (Vec<i32>, Vec<i32>) = edges.into_iter().unzip();
    let (xh, yh): (Vec<f64>, Vec<f64>) = holes.iter().map(|z| (z[0], z[1])).unzip();
    /*let tri = Builder::new()
        .set_points(xn, yn)
        .set_segments(xe, ye)
        .set_holes(xh, yh)
        .set_switches("Dnpqa.025")
        .build();*/
}
