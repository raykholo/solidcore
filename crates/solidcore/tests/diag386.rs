//! Diagnostic for #386 (not for merge): a torus as a surface of revolution of a
//! circle, its STEP domain adjusted to a v range starting below 0.
use solidcore::curve::curve::Curve;
use solidcore::curve::ellipse::{Ellipse, EllipseOptions};
use solidcore::geom::extent1d::Extent1d;
use solidcore::geom::extent2d::Extent2d;
use solidcore::geom::vector2d::Point2;
use solidcore::geom::vector3d::Vec3;
use solidcore::kernel::Kernel;
use solidcore::surface::surf_of_revolution::SurfOfRevolution;
use solidcore::surface::surface::Surface;

fn torus(k: &Kernel) -> SurfOfRevolution {
    let gen = Ellipse::new_circle(
        k,
        Vec3::new(10.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Extent1d::from_min_max(0.0, 360.0),
        2.0,
        &EllipseOptions::default(),
    )
    .unwrap();
    let mut s = SurfOfRevolution::new(
        k,
        Curve::from(gen),
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Extent2d::from_bounds(0.0, 360.0, 0.0, 360.0),
        false,
    )
    .unwrap();
    s.make_nurb(k).unwrap();
    s
}

#[test]
fn diag_386() {
    let k = Kernel::default();
    for (v0, v1) in [(-23.7, 156.3), (0.0, 180.0), (336.3, 516.3), (-360.0 - 23.7, -360.0 + 156.3)] {
        let mut s = torus(&k);
        let r = s.adjust_step_uv_domain(&k, &Extent2d::from_bounds(0.0, 360.0, v0, v1));
        let surf = Surface::from(s.clone());
        let d = surf.natural_uv_domain();
        let mut worst = 0.0f64;
        let n = 24;
        for i in 0..=n {
            for j in 0..=n {
                let u = d.min().x + (d.max().x - d.min().x) * f64::from(i) / f64::from(n);
                let v = d.min().y + (d.max().y - d.min().y) * f64::from(j) / f64::from(n);
                if let Ok(p) = surf.evaluate_point(Point2::new(u, v)) {
                    let rho = (p.x * p.x + p.y * p.y).sqrt();
                    let res = ((rho - 10.0).powi(2) + p.z * p.z).sqrt() - 2.0;
                    worst = worst.max(res.abs());
                }
            }
        }
        eprintln!(
            "DIAG386 v=[{v0},{v1}] adjust={:?} step={:?} natural={:?} worst_residual={worst:e}",
            r.as_ref().map(|_| ()),
            s.step_uv_domain(),
            d
        );
    }
}
