use ndarray::{Array1, Array2};
use crate::{parameters::Parameters, backend::{set_cavity_params,calc_dispersion}};
use gnuplot::*;

pub fn plot_disp_wrapper(mut prm:Parameters) {

    prm = set_cavity_params(prm);

    let omegas = Array1::linspace(prm.w_range.0, prm.w_range.1,prm.n_w);
    let q_pars = Array1::linspace(prm.q_range.0, prm.q_range.1, prm.n_q);

    let dispersion = calc_dispersion(&prm, &omegas, &q_pars);

    plot_disp(&prm,dispersion, omegas, q_pars);

}

fn plot_disp (prm: &Parameters,dispersion:Array2<f64>, omegas:Array1<f64>, q_pars: Array1<f64>) {
    let mut fig = Figure::new();

    let fname = format!("dispersion_q{}.png",prm.quality);

    fig.set_terminal("pngcairo size 1440,1080", fname.as_str())
    .set_pre_commands("set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb '#FFFFFF' behind");

    fig.axes2d()
    .set_x_range(AutoOption::Fix(prm.q_range.0), AutoOption::Fix(prm.q_range.1))
    .set_y_range(AutoOption::Fix(prm.w_range.0), AutoOption::Fix(prm.w_range.1))
    .set_x_ticks(Some((Auto,5)), &[MajorScale(2.0),OnAxis(false),Inward(false),Mirror(false)], &[Font("Helvetica", 24.0), TextColor("black")])
    .set_y_ticks(Some((Auto,5)), &[MajorScale(2.0),OnAxis(false),Inward(false),Mirror(false)], &[Font("Helvetica", 24.0), TextColor("black")])
    .set_cb_ticks(Some((Auto,5)),&[MajorScale(2.0),OnAxis(false),Inward(false),Mirror(false)], &[Font("Helvetica", 24.0), TextColor("black")])
    // .set_cb_label("", &[Font("Helvetica", 36.0), TextColor("white")])
    // .set_x_label("", &[Font("Helvetica", 36.0), TextColor("white")])
    // .set_y_label("", &[Font("Helvetica", 36.0), TextColor("white")])
    .set_margins(&[MarginLeft(0.10),MarginBottom(0.07), MarginRight(0.8)])
    .set_border(true, &[Bottom,Right,Left,Top], &[Color("black")])
    .image(dispersion, omegas.len(), q_pars.len(), Some((q_pars[0],omegas[0],q_pars.last().unwrap().clone(),omegas.last().unwrap().clone())), &[])
    ;

    let message = fig.save_to_png(fname, 1440, 1080);
    // let message_svg = fig.save_to_svg(fname, 1440, 1080);
    println!("{:?}", message);
    
}
