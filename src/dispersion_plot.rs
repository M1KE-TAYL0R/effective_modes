use std::f64::consts::PI;

use ndarray::{Array1, Array2};
use crate::{parameters::Parameters, backend::{set_cavity_params,calc_dispersion}};
use gnuplot::*;
use plotters::{prelude::*, style};
// use std::ops::Range;

pub fn plot_disp_wrapper(mut prm:Parameters) {

    prm = set_cavity_params(prm);
    
    prm.w_range = (prm.w_c * 0.95, prm.w_c * 1.2);


    let theta_max = 30.;

    let thetas = Array1::linspace(-theta_max, theta_max, prm.n_q);

    let q_pars = prm.w_c / prm.c * thetas.to_radians().tan();
    let _temp = q_pars.to_vec();
    let omegas = Array1::linspace(prm.w_range.0, prm.w_range.1,prm.n_w);
    // let q_pars = Array1::linspace(prm.q_range.0, prm.q_range.1, prm.n_q);

    let dispersion = calc_dispersion(&prm, &omegas, &q_pars);

    plot_disp(&prm,dispersion, omegas, thetas);
    // _plot_disp_plotters(&prm,dispersion, omegas, thetas).unwrap();

}

fn plot_disp (prm: &Parameters,dispersion:Array2<f64>, omegas:Array1<f64>, thetas: Array1<f64>) {

    // let omegas = omegas * 219474.63;

    let log_min = 1.0;
    let log_max = 6.0;

    let omegas = omegas / prm.w_c;

    let mut fig = Figure::new();

    let fname = format!("dispersion_q{}.png",prm.quality);

    fig.set_terminal("pngcairo size 1250,1080", fname.as_str())
    .set_pre_commands("set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb '#FFFFFF' behind");

    let log_disp = dispersion.log10();

    fig.axes2d()
    .set_x_range(AutoOption::Fix(thetas[0]), AutoOption::Fix(thetas.last().unwrap().clone()))
    .set_y_range(AutoOption::Fix(omegas[0]), AutoOption::Fix(omegas.last().unwrap().clone()))
    .set_x_ticks(Some((AutoOption::Fix(20.),0)), &[MajorScale(4.0),OnAxis(false),Inward(false),Mirror(false),MinorScale(2.0)], &[Font("Helvetica", 48.0), TextColor("black"), TextOffset(-0.85,-1.0)])
    .set_y_ticks(Some((AutoOption::Fix(0.1),5)), &[MajorScale(4.0),OnAxis(false),Inward(false),Mirror(false),MinorScale(2.0)], &[Font("Helvetica", 48.0), TextColor("black")])
    .set_cb_ticks(Some((Auto,5)),&[MajorScale(5.0),OnAxis(false),Inward(false),Mirror(false),MinorScale(2.0)], &[Font("Helvetica", 24.0), TextColor("black")])
    .set_cb_range(AutoOption::Fix(log_min), AutoOption::Fix(log_max))
    // .set_cb_label("", &[Font("Helvetica", 36.0), TextColor("white")])
    // .set_x_label("", &[Font("Helvetica", 36.0), TextColor("white")])
    // .set_y_label("", &[Font("Helvetica", 36.0), TextColor("white")])
    .set_margins(&[MarginLeft(0.10),MarginBottom(0.11), MarginRight(0.85), MarginTop(0.95)])
    .set_border(true, &[Bottom,Right,Left,Top], &[Color("black"),LineWidth(3.0)])
    .image(log_disp, omegas.len(), thetas.len(), Some((thetas[0],omegas[0],thetas.last().unwrap().clone(),omegas.last().unwrap().clone())), &[])
    ;

    let message = fig.save_to_png(fname, 1440, 1080);
    // let message_svg = fig.save_to_svg(fname, 1440, 1080);
    println!("{:?}", message);
    
}

fn _plot_disp_plotters(prm: &Parameters,dispersion:Array2<f64>, omegas:Array1<f64>, q_pars: Array1<f64>) -> Result<(), Box<dyn std::error::Error>>{
    
    let fname = format!("dispersion_q{}.png",prm.quality);

    let root = BitMapBackend::new(fname.as_str(), (1080, 1080)).into_drawing_area();
    root.fill(&WHITE)?;

    let x_spec = q_pars[0]..q_pars.last().unwrap().clone();
    let y_spec = omegas[0]..omegas.last().unwrap().clone();
    let dq = (prm.q_range.1 - prm.q_range.0) / (prm.n_q as f64);
    let dw = (prm.w_range.1 - prm.w_range.0) / (prm.n_w as f64);

    let original_style = ShapeStyle {
        color: BLACK.into(),
        filled: false,
        stroke_width: 4,
    };
    
    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .x_label_area_size(10)
        .y_label_area_size(10)
        .build_cartesian_2d(x_spec, y_spec)?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .set_all_tick_mark_size(30)
        .axis_style(original_style)
        .draw()?;

    // let plotting_area = chart.plotting_area();

    // let range = plotting_area.get_pixel_range();

    // let (pw, ph) = (range.0.end - range.0.start, range.1.end - range.1.start);
    // let (xr, yr) = (chart.x_range(), chart.y_range());

    // dispersion.indexed_iter().map(|((q_ind,w_ind),val)| {
    //     let w = omegas[w_ind];
    //     let q = q_pars[q_ind];
    // });

    let _max_val = *dispersion.iter().max_by(|&x,&y| x.total_cmp(y)).unwrap();

    let mut min_val = 100000.;

    chart.draw_series(
        dispersion
            .indexed_iter()
            .map(|((w_ind,q_ind),val)| {
                let w = omegas[w_ind];
                let q = q_pars[q_ind];

                if *val < min_val {
                    min_val = *val;
                }

                Rectangle::new(
                    [(q, w), (q + dq, w + dw)],
                    // VulcanoHSL::get_color_normalized(*val as f32, 0.0 as f32, max_val as f32),
                    _get_pm3d_color_norm(val.log10(), 6., -2.0)
                )
            }),
    )?;

    println!("Min R = {}", min_val);
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

fn _get_pm3d_color_norm(val:f64, max:f64, min:f64) -> style::RGBColor {
    
    let norm_val = (val - min) / (max - min);

    let r = (255.0 * norm_val.sqrt()).floor() as u8;
    let g = (255.0 * norm_val.powi(3)).floor() as u8;

    let mut b = 0 as u8;

    if (2. * PI * norm_val).sin() > 0.0 {
        b = (255.0 * (2. * PI * norm_val).sin()) as u8;
    }

    RGBColor(r, g, b)
}