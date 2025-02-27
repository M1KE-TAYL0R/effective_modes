use std::f64::consts::PI;
use ndarray::{Array1, Array2, s};
use gnuplot::*;
use plotters::prelude::*;
use statrs::statistics::Statistics;
use ndarray_linalg::c64;

fn main() {
    let mut prm = Parameters{
        l_c: 0.0,
        c: 1.0 / 137.0,
        r: 0.95,
        w_0: 1.0,
        n_w: 5000,
        n_q: 1000,
        n_w_bins: 500,
        del_k: 1.0,
        quality: 5.0,
        q_range: (0.0,100.0),
        w_range: (0.95, 1.05)
    };

    // let gamma = 4.0 * prm.w_0 * (1.0 - prm.r) / prm.r.sqrt();

    if prm.quality == 0.0 {
        let gamma = - prm.w_0 / (2.0 * PI) * prm.r.powi(4).ln() ;
        prm.del_k = gamma / prm.c;
        prm.quality = prm.w_0 / gamma;
    }
    else {
        let gamma = prm.w_0 / prm.quality;
        prm.del_k = gamma / prm.c;
        prm.r = (- 2.0 * PI / prm.quality).exp().powf(0.25);
        println!("Reflectance: {}", prm.r);
    }


    println!("Quality Factor: {}", prm.quality);
    println!("Delta q_perp: {}", prm.del_k);


    let q_0 = prm.w_0/prm.c;
    prm.l_c  = 2.0 * PI / q_0;

    let omegas = Array1::linspace(prm.w_range.0, prm.w_range.1,prm.n_w);
    let q_pars = Array1::linspace(prm.q_range.0, prm.q_range.1, prm.n_q);

    let mut disp: Array2<f64> = Array2::zeros((prm.n_q,prm.n_w));

    disp.indexed_iter_mut().for_each(|(ind, val)| {
        let q_par = q_pars[ind.0];
        let omega = omegas[ind.1];

        *val = enhancement_function(omega, q_par, &prm);
    });

    let dispersion = disp.t().to_owned();

    let weights = bin_dispersion(&prm, omegas.clone(), q_pars.clone(), dispersion.clone());
    plot_weights(weights, &prm, &omegas).unwrap();

    // println!("test: {}", enhancement_function(1.0, 0.0, &prm));

    // plot_disp(dispersion, omegas, q_pars);
}

fn plot_weights(weights:Array1<f64>, prm:&Parameters, omegas:&Array1<f64>) -> Result<(), Box<dyn std::error::Error>> {
    let omegas = omegas.slice(s![..;prm.n_w / prm.n_w_bins]);

    let x = omegas.to_vec();
    let y: Vec<f64> = weights.to_vec();
    let data: Vec<(f64,f64)> = x.into_iter().zip(y).collect();

    let root = SVGBackend::new("weights.svg", (1440,1080)).into_drawing_area();
        
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        // .caption("Test Dispersion", ("helvetica", 50*scale_factor))
        .x_label_area_size(70)
        .y_label_area_size(100)
        .build_cartesian_2d(prm.w_range.0..prm.w_range.1, 0.0 .. weights.max())?;

    chart.configure_mesh()
    .x_label_style(("helvetica", 50))
    .y_label_style(("helvetica", 50))
    // .disable_mesh()
    .set_all_tick_mark_size(20)
    // .bold_line_style(original_style)
    // .x_label_formatter(&|v| format!("{:.1}", v))
    // .y_label_formatter(&|v| format!("{:.1}", v))
    .draw()?;

    chart.draw_series(LineSeries::new((data.clone()).into_iter().map(|(x,y)| (x, y)), BLUE))?;

    // println!("{:?}",data);
    
    Ok(())
}

fn integrate_bin (mut bin:Array2<f64>, omegas_bin:Array1<f64>, q_pars: &Array1<f64>, prm:&Parameters) -> f64 {
    
    let dq = (prm.q_range.1 - prm.q_range.0) / (prm.n_q as f64);
    let dw = (prm.w_range.1 - prm.w_range.0) / (prm.n_w as f64);

    bin.indexed_iter_mut().for_each(|(ind, val)| {
        let q_par = q_pars[ind.1];
        let omega = omegas_bin[ind.0];
        let q_z = (omega.powi(2) - prm.c*prm.c * q_par.powi(2)).sqrt();
        
        // let jacobian = 1.0; 
        let mut jacobian =  q_par * omega / prm.c / q_z; //.powf(-1.0);

        
        let qn_2 = (omega / prm.c).powf(2.0);
        let qn_z_2 = (2.0 * PI / prm.l_c).powf(2.0);
        let qn_par = (qn_2 - qn_z_2).sqrt();
        // let qn_par = q_par;
        
        if prm.del_k < qn_par {

            let del_theta = 2.0 * (prm.del_k / qn_par).asin();
            // let round_del_t = ( PI / del_theta).ceil();

            // jacobian *= PI / round_del_t;=

            // if del_theta > 2. * PI {
            //     del_theta = 2. * PI;
            // }
            jacobian *= del_theta;
        }
        else {
            jacobian *= PI;
        }

        // if (*val * jacobian).is_finite() {
        *val = val.clone() * jacobian;
        // }
        // else {
        //     println!("Omega = {}, q_par = {}, q_z = {}",omega,q_par,q_z);
        //     *val = 0.0;
        // }
    });

   let weight = bin.sum() * dq* dw;

   weight
}

fn bin_dispersion (prm:&Parameters, omegas:Array1<f64>, q_pars: Array1<f64>, dispersion:Array2<f64>) -> Array1<f64>{
    
    let bin_size = prm.n_w / prm.n_w_bins;
    let mut weights: Array1<f64> = Array1::zeros( prm.n_w_bins);

    for bin_ind in 0..prm.n_w_bins {
        let omegas_bin = omegas.slice(s![(bin_ind * bin_size) .. ((bin_ind+1)*bin_size)]);
        let bin = dispersion.slice(s![(bin_ind * bin_size) .. ((bin_ind+1)*bin_size),..]);

        // println!("{:?}", omegas_bin.shape());

        weights[bin_ind] = integrate_bin(bin.to_owned(), omegas_bin.to_owned(), &q_pars, prm);
    }
    
    weights
}

fn _plot_disp (dispersion:Array2<f64>, omegas:Array1<f64>, q_pars: Array1<f64>) {
    let mut fig = Figure::new();

    let fname = "test.png";

    fig.set_terminal("pngcairo size 1440,1080", fname)
    .set_pre_commands("set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb '#000000' behind");

    fig.axes2d()
    // .set_x_range(AutoOption::Fix(x_min), AutoOption::Fix(-x_min))
    // .set_y_range(AutoOption::Fix(0.0), AutoOption::Fix(prm.max_energy))
    .set_x_ticks(Some((Auto,5)), &[MajorScale(2.0),OnAxis(false),Inward(false),Mirror(false)], &[Font("Helvetica", 24.0), TextColor("white")])
    .set_y_ticks(Some((Auto,5)), &[MajorScale(2.0),OnAxis(false),Inward(false),Mirror(false)], &[Font("Helvetica", 24.0), TextColor("white")])
    .set_cb_ticks(Some((Auto,5)),&[MajorScale(2.0),OnAxis(false),Inward(false),Mirror(false)], &[Font("Helvetica", 24.0), TextColor("white")])
    // .set_cb_label("", &[Font("Helvetica", 36.0), TextColor("white")])
    // .set_x_label("", &[Font("Helvetica", 36.0), TextColor("white")])
    // .set_y_label("", &[Font("Helvetica", 36.0), TextColor("white")])
    // .set_margins(&[MarginLeft(0.07),MarginBottom(0.07), MarginRight(0.87)])
    // .set_border(true, &[Bottom,Right,Left,Top], &[Color("white")])
    .image(dispersion, omegas.len(), q_pars.len(), Some((q_pars[0],omegas[0],q_pars.last().unwrap().clone(),omegas.last().unwrap().clone())), &[])
    ;

    let message = fig.save_to_png(fname, 1440, 1080);
    // let message_svg = fig.save_to_svg(fname, 1440, 1080);
    println!("{:?}", message);
    
}


fn enhancement_function (omega:f64, q_par:f64, prm:&Parameters) -> f64 {

    let q_z = (omega.powi(2) / prm.c.powi(2) - q_par.powi(2)).sqrt();

    if q_z <  PI / prm.l_c || q_z.is_nan(){
        return 0.0;
    }
    // let jacobian = omega / prm.c.powi(2);
    // let t = 1.0 - prm.r;
    let r = prm.r;

    let mut weight: c64 = 1.0 + r * (prm.l_c * q_z * c64::i()).exp(); 
    // let mut weight: c64 = t + t * r * (prm.l_c * q_z * c64::i()).exp(); 
    weight /= 1.0 - r*r*(2.0 * prm.l_c * q_z * c64::i()).exp();
    // weight -= 1.0;

    // let mut weight = t.powi(2) / (1.0 + prm.r.powi(2) - 2.0 * prm.r * (prm.l_c * q_z).cos());
    // weight += 1.0;
    // weight += t * 2. * (r.powi(3) * (3. * prm.l_c * q_z).cos() + r.powi(2) * (2. * prm.l_c * q_z).cos() - r * (prm.l_c * q_z).cos() - 1.) 
    //             / ( r.powi(4) -  2. * r.powi(2) * (2. * prm.l_c * q_z).cos() + 1.);

    let re_weight = weight.norm() - 1.0;

    if re_weight < 0.0 { return 0.0}

    re_weight * re_weight
}

pub struct Parameters {
    pub l_c: f64,
    pub c: f64,
    pub r: f64,
    pub w_0: f64,
    pub n_w: usize,
    pub n_q: usize,
    pub n_w_bins: usize,
    pub del_k: f64,
    pub quality: f64,
    pub q_range: (f64,f64),
    pub w_range: (f64,f64)
}