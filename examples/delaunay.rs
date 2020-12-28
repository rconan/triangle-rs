use triangle_rs as dtri;
use triangle_rs::TriDraw;

fn main() {
    let p_x: Vec<f64> = vec![0., 1., 0., -1., 0., 1.];
    let p_y: Vec<f64> = vec![0., 0., 1., 0., -1.5, 1.5];
    let tri = dtri::Builder::new().set_points(p_x, p_y).build();
    println!("Delaunay: {:#?}", tri);

    tri.mesh("examples/delaunay.svg", 2f64);
    /*
    let plot = BitMapBackend::new("examples/delaunay.png", (768, 768)).into_drawing_area();
    plot.fill(&WHITE).unwrap();
    let mut chart = ChartBuilder::on(&plot)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(-2.0..2.0, -2.0..2.0)
        .unwrap();
    chart.configure_mesh().draw().unwrap();
    chart
        .draw_series(
            p_x.iter()
                .cloned()
                .zip(p_y.iter().cloned())
                .map(|p| Circle::new(p, 5, &RED)),
        )
        .unwrap();
    let vertices: Vec<Vec<(f64, f64)>> = tri
        .triangles
        .iter()
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
    */
}
