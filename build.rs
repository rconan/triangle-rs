use cc;

fn main() {
    cc::Build::new()
        .file("src/triangle.c")
        .flag("-O")
        .flag("-DLINUX")
        .flag("-DTRILIBRARY")
        .include("src")
        .define("REAL", "double")
        .compile("triangle");
}
