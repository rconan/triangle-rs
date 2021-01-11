use complot as plt;
use triangle_rs::Builder;
use plt::TriPlot;

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
    tri.mesh(&tri.x(), &tri.y(), [0; 3], &mut ax);
}
