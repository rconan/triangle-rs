use complot as plt;
use complot::TriPlot;
use spade::delaunay::FloatDelaunayTriangulation;
use triangle_rs as dtri;
/*fn fun(x: f64, y: f64) -> f64 {
    y.powf(3.) - 2. * x.powf(2.)
}*/

fn main() {
    let nodes = vec![1.5,1.5,1.,0.,0.,1.,-1.,0.,0.,-1.5,0.,0.,-0.5,-0.5];
    let mut builder = dtri::Builder::new();
    let tri =builder.add_nodes(&nodes).build();

    println!("Delaunay: {:#?}", tri);

    /*    (0..tri.triangles.len())
            .for_each(|k| println!("TRI#{}: {}", k, tri.winding_number(&[0., 0.], k)));
    */
    println!(
        "Found in triangle #{:?}",
        tri.which_contains_point(&[0.6, 0.5])
    );

    let fig = plt::canvas("examples/delaunay.svg");
    let mut ax = plt::chart([-2., 2., -2., 2.], &fig);
    let (p_x,p_y): (Vec<f64>,Vec<f64>) = nodes.chunks(2).map(|xy| (xy[0],xy[1])).unzip();
    plt::trimesh(&p_x, &p_y, [0; 3], &mut ax);

    /*
    let z: Vec<f64> = p_x
        .iter()
        .zip(p_y.iter())
        .map(|(x, y)| fun(*x, *y))
        .collect();
    let p = vec![-0.5, -0.5];
    let zi = tri.barycentric_interpolation(&p, &z);
    println!("{:.6}/{:.6}", fun(p[0], p[1]), zi);

    let q = vec![&[0., 0.], &[1., 1.]];
    let qf: Vec<f64> = q.into_iter().flatten().cloned().collect();
    */
}
