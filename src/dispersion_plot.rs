use ndarray::{Array1, Array2};
use crate::{parameters::Parameters, backend::{set_cavity_params,calc_dispersion}};
use gnuplot::*;
use plotters::prelude::*;
// use std::ops::Range;

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
    .set_x_ticks(Some((Auto,5)), &[MajorScale(5.0),OnAxis(false),Inward(false),Mirror(false)], &[Font("Helvetica", 24.0), TextColor("black")])
    .set_y_ticks(Some((Auto,5)), &[MajorScale(5.0),OnAxis(false),Inward(false),Mirror(false)], &[Font("Helvetica", 24.0), TextColor("black")])
    .set_cb_ticks(Some((Auto,5)),&[MajorScale(5.0),OnAxis(false),Inward(false),Mirror(false)], &[Font("Helvetica", 24.0), TextColor("black")])
    // .set_cb_label("", &[Font("Helvetica", 36.0), TextColor("white")])
    // .set_x_label("", &[Font("Helvetica", 36.0), TextColor("white")])
    // .set_y_label("", &[Font("Helvetica", 36.0), TextColor("white")])
    .set_margins(&[MarginLeft(0.10),MarginBottom(0.07), MarginRight(0.8)])
    .set_border(true, &[Bottom,Right,Left,Top], &[Color("black"),LineWidth(3.0)])
    .image(dispersion, omegas.len(), q_pars.len(), Some((q_pars[0],omegas[0],q_pars.last().unwrap().clone(),omegas.last().unwrap().clone())), &[])
    ;

    let message = fig.save_to_png(fname, 1440, 1080);
    // let message_svg = fig.save_to_svg(fname, 1440, 1080);
    println!("{:?}", message);
    
}

fn _plot_disp_plotters(prm: &Parameters,dispersion:Array2<f64>, omegas:Array1<f64>, q_pars: Array1<f64>) -> Result<(), Box<dyn std::error::Error>>{
    
    let fname = format!("dispersion_q{}.png",prm.quality);

    let root = BitMapBackend::new(fname.as_str(), (1440, 1080)).into_drawing_area();
    root.fill(&WHITE)?;

    let x_spec = q_pars[0]..q_pars.last().unwrap().clone();
    let y_spec = omegas[0]..omegas.last().unwrap().clone();
    let dq = (prm.q_range.1 - prm.q_range.0) / (prm.n_q as f64);
    let dw = (prm.w_range.1 - prm.w_range.0) / (prm.n_w as f64);

    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .x_label_area_size(10)
        .y_label_area_size(10)
        .build_cartesian_2d(x_spec, y_spec)?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .draw()?;

    // let plotting_area = chart.plotting_area();

    // let range = plotting_area.get_pixel_range();

    // let (pw, ph) = (range.0.end - range.0.start, range.1.end - range.1.start);
    // let (xr, yr) = (chart.x_range(), chart.y_range());

    // dispersion.indexed_iter().map(|((q_ind,w_ind),val)| {
    //     let w = omegas[w_ind];
    //     let q = q_pars[q_ind];
    // });

    let max_val = *dispersion.iter().max_by(|&x,&y| x.total_cmp(y)).unwrap();

    chart.draw_series(
        dispersion
            .indexed_iter()
            .map(|((q_ind,w_ind),val)| {
                let w = omegas[w_ind];
                let q = q_pars[q_ind];

                Rectangle::new(
                    [(q, w), (q + dq, w + dw)],
                    VulcanoHSL::get_color_normalized(*val as f32, 0.0 as f32, max_val as f32),
                )
            }),
    )?;

    // let mut matrix = [[0; 15]; 15];

    // for i in 0..15 {
    //     matrix[i][i] = i + 4;
    // }

    // chart.draw_series(
    //     matrix
    //         .iter()
    //         .zip(0..)
    //         .flat_map(|(l, y)| l.iter().zip(0..).map(move |(v, x)| (x, y, v)))
    //         .map(|(x, y, v)| {
    //             Rectangle::new(
    //                 [(x, y), (x + 1, y + 1)],
    //                 HSLColor(
    //                     240.0 / 360.0 - 240.0 / 360.0 * (*v as f64 / 20.0),
    //                     0.7,
    //                     0.1 + 0.4 * *v as f64 / 20.0,
    //                 )
    //                 .filled(),
    //             )
    //         }),
    // )?;


    Ok(())
}