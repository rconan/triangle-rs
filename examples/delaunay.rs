use triangle_rs as dtri;
use triangle_rs::TriDraw;

fn fun(x: f64, y: f64) -> f64 {
    y.powf(3.) - 2.*x.powf(2.)
}

fn main() {
    let p_x: Vec<f64> = vec![1.5, 1., 0., -1.,  0. , 0.];
    let p_y: Vec<f64> = vec![1.5, 0., 1.,  0., -1.5, 0.];
    let tri = dtri::Builder::new().set_points(p_x.clone(), p_y.clone()).build();
    println!("Delaunay: {:#?}", tri);

/*    (0..tri.triangles.len())
        .for_each(|k| println!("TRI#{}: {}", k, tri.winding_number(&[0., 0.], k)));
*/
    println!("Found in triangle #{:?}",tri.which_contains_point(&[0.6,0.5]));

    tri.mesh("examples/delaunay.svg", 2f64);

    let z: Vec<f64> = p_x.iter().zip(p_y.iter()).map(|(x,y)| fun(*x,*y)).collect();
    let p = vec![-0.5,-0.5];
    let zi = tri.barycentric_interpolation(&p, &z);
    println!("{:.6}/{:.6}",fun(p[0],p[1]),zi);

    let q = vec![&[0.,0.],&[1.,1.]];
    let qf: Vec<f64> = q.into_iter().flatten().cloned().collect();
}
