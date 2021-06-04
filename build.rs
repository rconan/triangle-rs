use cc;

fn main() {
    println!("cargo:rerun-if-changed=src/triangle.c");
    cc::Build::new()
        .file("src/triangle.c")
        .flag("-O")
        .flag("-DLINUX")
        .flag("-DTRILIBRARY")
        .include("src")
        .define("REAL", "double")
        .warnings(false)
        .compile("triangle");
}
