use std::f64::consts::PI;
use ndarray::{Array1, Array2};
use gnuplot::*;

fn main() {
    let mut prm = Parameters{
        l_c: 0.0,
        c: 1.0 / 137.0,
        r: 0.9,
        n_w_bins: 1000,
        n_q_bins: 1000
    };

    let w_0 = 1.0;
    let q_0 = w_0/prm.c;
    prm.l_c  = 2.0 * PI / q_0;

    let omegas = Array1::linspace(1.0, 1.5,prm.n_w_bins);
    let q_pars = Array1::linspace(0.0, 1.0, prm.n_q_bins);

    let mut dispersion: Array2<f64> = Array2::zeros((prm.n_w_bins,prm.n_q_bins));

    dispersion.indexed_iter_mut().for_each(|(ind, val)| {
        let q_par = q_pars[ind.1];
        let omega = omegas[ind.0];

        *val = enhancement_function(omega, q_par, &prm);
    });

    println!("test: {}", enhancement_function(1.0, 0.0, &prm));

    // _plot_disp(dispersion, omegas, q_pars);
}

fn _plot_disp (dispersion:Array2<f64>, omegas:Array1<f64>, q_pars: Array1<f64>) {
    let mut fig = Figure::new();

    let fname = "test.png";

    fig.set_terminal("pngcairo size 1440,1080", fname)
    .set_pre_commands("set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb '#000000' behind");

    fig.axes2d()
    // .set_x_range(AutoOption::Fix(x_min), AutoOption::Fix(-x_min))
    // .set_y_range(AutoOption::Fix(0.0), AutoOption::Fix(prm.max_energy))
    // .set_x_ticks(Some((Auto,5)), &[MajorScale(2.0),OnAxis(false),Inward(false),Mirror(false)], &[Font("Helvetica", 24.0), TextColor("white")])
    // .set_y_ticks(Some((Auto,5)), &[MajorScale(2.0),OnAxis(false),Inward(false),Mirror(false)], &[Font("Helvetica", 24.0), TextColor("white")])
    // .set_cb_ticks(Some((Auto,5)),&[MajorScale(2.0),OnAxis(false),Inward(false),Mirror(false)], &[Font("Helvetica", 24.0), TextColor("white")])
    // .set_cb_label("", &[Font("Helvetica", 36.0), TextColor("white")])
    // .set_x_label("", &[Font("Helvetica", 36.0), TextColor("white")])
    // .set_y_label("", &[Font("Helvetica", 36.0), TextColor("white")])
    // .set_margins(&[MarginLeft(0.07),MarginBottom(0.07), MarginRight(0.87)])
    // .set_border(true, &[Bottom,Right,Left,Top], &[Color("white")])
    .image(dispersion, q_pars.len(), omegas.len(), Some((0.0,0.0,q_pars.last().unwrap().clone(),omegas.last().unwrap().clone())), &[])
    ;

    let message = fig.save_to_png(fname, 1440, 1080);
    // let message_svg = fig.save_to_svg(fname, 1440, 1080);
    println!("{:?}", message);
    
}


fn enhancement_function (omega:f64, q_par:f64, prm:&Parameters) -> f64 {

    let q_z = (omega.powi(2) / prm.c.powi(2) + q_par.powi(2)).sqrt();
    // let jacobian = omega / prm.c.powi(2);
    let t = 1.0 - prm.r;

    t.powi(2) / (1.0 + prm.r.powi(2) - 2.0 * prm.r * (prm.l_c * q_z).cos())
}

pub struct Parameters {
    pub l_c: f64,
    pub c: f64,
    pub r: f64,
    pub n_w_bins: usize,
    pub n_q_bins: usize
}