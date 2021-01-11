use complot as plt;
use plotters::prelude::*;
use triangle_rs::Builder;

fn main() {
    let p0 = [1., 1., -1., 1., -1., -1., 1., -1.];
    let p1 = [0.5, 0., 0., 0.5, -0.5, 0., 0., -0.5];
    let mut builder = Builder::new();
    builder
        .add_polygon(&p0)
        .add_polygon(&p1)
        .add_holes(0., 0.)
        .set_switches("pDqa0.01");
    let tri = builder.build();
    let fig = plt::canvas("examples/box.svg");
    let mut ax = plt::chart([-2., 2., -2., 2.], &fig);
    tri.triangle_iter()
        .map(|t| {
            t.iter()
                .map(|&i| (tri.x()[i], tri.y()[i]))
                .collect::<Vec<(f64, f64)>>()
        })
        .for_each(|v| {
            ax.draw_series(LineSeries::new(
                v.iter().cycle().take(4).map(|(x, y)| (*x, *y)),
                &BLACK,
            ))
            .unwrap();
        });
}
