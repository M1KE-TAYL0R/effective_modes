use std::f64::consts::PI;
use ndarray::{s, Array1, Array2};
use gnuplot::*;
use plotters::{prelude::*, style::{full_palette::{DEEPORANGE_A400, LIGHTBLUE_A700, LIME_A700}, Color}};
use statrs::statistics::Statistics;
use ndarray_linalg::c64;
use ndarray_npy::{write_npy,read_npy};
use std::path::Path;
use rayon::prelude::*;


fn define_parameters() -> Parameters {
    let prm = Parameters{
        l_c: 0.0,                         // Placeholder
        c: 1.0 / 137.0,                   // Speed of light in atomic units
        r: 0.0,                           // Placeholder if quality != 0
        w_0: 0.1,                         // Cavity resonance frequency in atomic units
        n_w: 1000,                       // Number of \omega grid points
        n_q: 10000,                      // Number of q_\| grid points per bin integration
        n_w_bins: 5000,                   // Number of omega_n bins
        del_k: 1.0,                       // Value of \Delta q_\perp
        quality: 500.0,                   // Cavity Quality Factor
        q_range: (0.0,10.0),              // Range of q_\| points integrated over
        w_range: (0.09, 0.12),            // Range of omega_n
        routine: "ARDOF".to_string()      // Set the routine: ManyQ, SingQ, Dispn, ARDOF
    };

    prm
}


fn main() {
    // Set Parameters
    let prm = define_parameters();

    // Run routine specified by prm.routine
    match prm.routine.as_str() {
        "ManyQ"  => many_q_factors(prm),
        "SingQ"  => single_q_factor(prm),
        "Dispn"  => plot_disp_wrapper(prm),
        "ARDOF"  => ardof(prm),
        _        => panic!("Invalid Routine Parameter: check prm.routine")
    };
}

fn ardof(mut prm: Parameters) {
    let q_0 = prm.w_0/prm.c;
    prm.l_c  = 2.0 * PI / q_0;

    let omegas = Array1::linspace(prm.w_range.0, prm.w_range.1,prm.n_w);
    let q_pars = Array1::linspace(prm.q_range.0, prm.q_range.1, prm.n_q);

    let qualities = vec![5000.,500.,50.];

    let mut ardofs:Vec<(f64, Array1<f64>)> = Vec::new();

    let dq = (prm.q_range.1 - prm.q_range.0) / (prm.n_q as f64);

    qualities.iter().enumerate().for_each(|(ind, quality)| {

        prm.quality = *quality;

        ardofs.push((prm.quality,Array1::zeros(prm.n_w)));

        println!("Calculating Q = {}", quality);

        // Define r
        let gamma = prm.w_0 / prm.quality;
        prm.del_k = gamma / prm.c;

        // Lorentzian linewidth approx
        prm.r = (- 2.0 * PI / prm.quality).exp().powf(0.25);

        // let mut ardof: Array1<f64> = Array1::zeros([prm.n_w]);

        omegas.iter().enumerate().for_each(|(ind1, omega)|  {

            // let mut val = 0.0;
            // q_pars.iter().for_each(|q_par| {
            let mut temp: Vec<f64> = Array1::zeros(prm.n_q).to_vec();
            q_pars.to_vec().into_par_iter().zip(temp.par_iter_mut()).for_each(|(q_par,temp1)| {

                let q_z = (omega.powi(2) - prm.c*prm.c * q_par.powi(2)).sqrt();
                
                let jacobian =  q_par * omega / prm.c / q_z;

                *temp1 = enhancement_function(*omega, q_par, &prm) * jacobian;
            });
            ardofs[ind].1[ind1] =  temp.iter().sum();
            ardofs[ind].1[ind1] *= dq;

        });

        // ardofs[ind].1 *= dq;
    });

    plot_ardofs_quality(&ardofs, &prm, &omegas).unwrap();
}

fn set_cavity_params(mut prm:Parameters) -> Parameters {
    // Calculate prm.del_k and prm.r from prm.quality
    if prm.quality == 0.0 {
        // Lorentzian linewidth approx
        let gamma = - prm.w_0 / (2.0 * PI) * prm.r.powi(4).ln(); 

        prm.del_k = gamma / prm.c;
        prm.quality = prm.w_0 / gamma;
    }
    // Calculate prm.del_k and prm.quality from prm.r 
    else {
        let gamma = prm.w_0 / prm.quality;
        prm.del_k = gamma / prm.c;

        // Lorentzian linewidth approx
        prm.r = (- 2.0 * PI / prm.quality).exp().powf(0.25);
    }

    println!("Reflectance: {}", prm.r);
    println!("Quality Factor: {}", prm.quality);
    println!("Delta q_perp: {}", prm.del_k);

    let q_0 = prm.w_0/prm.c;
    prm.l_c  = 2.0 * PI / q_0;

    prm
}

fn plot_disp_wrapper(mut prm:Parameters) {

    prm = set_cavity_params(prm);

    let omegas = Array1::linspace(prm.w_range.0, prm.w_range.1,prm.n_w);
    let q_pars = Array1::linspace(prm.q_range.0, prm.q_range.1, prm.n_q);

    let dispersion = calc_dispersion(&prm, &omegas, &q_pars);

    plot_disp(dispersion, omegas, q_pars);

}

fn many_q_factors(mut prm:Parameters) {
    let q_0 = prm.w_0/prm.c;
    prm.l_c  = 2.0 * PI / q_0;

    let omegas = Array1::linspace(prm.w_range.0, prm.w_range.1,prm.n_w);
    let q_pars = Array1::linspace(prm.q_range.0, prm.q_range.1, prm.n_q);

    let qualities = vec![5000.,500.,50.];

    let (prm,output) = scan_quality_factors(qualities, prm, omegas, q_pars);

    save_output(&prm,output);
}

fn save_output (prm: &Parameters, output: Vec<(f64,Array1<f64>)>) {

    output.iter().for_each(|(qual,weight)| {
        let fname = format!("ell_n__w{}_nq{}_qual{}.npy", prm.w_0, prm.n_q, qual);
        println!("Saving: '{}'", fname);

        write_npy(fname, weight).unwrap();
    });
}

fn single_q_factor(mut prm:Parameters) {

    prm = set_cavity_params(prm);

    let omegas = Array1::linspace(prm.w_range.0, prm.w_range.1,prm.n_w);
    let q_pars = Array1::linspace(prm.q_range.0, prm.q_range.1, prm.n_q);

    let dispersion = calc_dispersion(&prm, &omegas, &q_pars);

    let weights = bin_dispersion(&prm, omegas.clone(), q_pars.clone(), dispersion.clone());
    plot_weights(weights, &prm, &omegas).unwrap();

}

fn scan_quality_factors(qualities:Vec<f64>, mut prm:Parameters, omegas:Array1<f64>, q_pars: Array1<f64>) -> (Parameters, Vec<(f64,Array1<f64>)>){

    let mut weights_quality:Vec<Array1<f64>> = Vec::new();

    qualities.iter().for_each(|quality| {
        prm.quality = *quality;

        println!("Calculating Q = {}", quality);

        let fname = format!("ell_n__w{}_nq{}_qual{}.npy", prm.w_0, prm.n_q, prm.quality);

        if Path::new(fname.as_str()).is_file() {
            let weights: Array1<f64> = read_npy(fname).unwrap();
            weights_quality.push(weights);
        }
        else {
             // Define r
             let gamma = prm.w_0 / prm.quality;
             prm.del_k = gamma / prm.c;
 
             // Lorentzian linewidth approx
             prm.r = (- 2.0 * PI / prm.quality).exp().powf(0.25);

            let dispersion = calc_dispersion(&prm, &omegas, &q_pars);
 
            weights_quality.push(bin_dispersion(&prm, omegas.clone(), q_pars.clone(), dispersion));
        }
    });

    plot_weights_quality(&weights_quality, &prm, &omegas, &qualities).unwrap();

    (prm, qualities.into_iter().zip(weights_quality).collect())

}

fn calc_dispersion (prm:&Parameters, omegas:&Array1<f64>, q_pars: &Array1<f64>) -> Array2<f64>{
    let mut disp: Array2<f64> = Array2::zeros((prm.n_q,prm.n_w));

    disp.indexed_iter_mut().for_each(|(ind, val)| {
        let q_par = q_pars[ind.0];
        let omega = omegas[ind.1];

        *val = enhancement_function(omega, q_par, &prm);
    });

    disp.t().to_owned()
}

fn plot_ardofs_quality(ardofs:&Vec<(f64,Array1<f64>)>, prm: &Parameters, omegas:&Array1<f64>) ->  Result<(), Box<dyn std::error::Error>> {
    
    // Find maximum y-value
    let max_q_ind = ardofs.iter()
        .enumerate()
        .max_by(|(_, (a,_)), (_, (b, _))| a.total_cmp(b))
        .map(|(index, _)| index);
    // let test = ardofs[max_q_ind.unwrap()].to_vec();
    let max_y = *ardofs[max_q_ind.unwrap()].1
        .to_vec()
        .iter()
        .max_by(|a,b| a.total_cmp(b))
        .unwrap();

    println!("Max y = {}", max_y);

    let min_y = *ardofs[max_q_ind.unwrap()].1
        .to_vec()
        .iter()
        .min_by(|a,b| a.total_cmp(b))
        .unwrap();

    let root = SVGBackend::new("ardofs.svg", (1440,1080)).into_drawing_area();
    
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(40)
        // .caption("Test Dispersion", ("helvetica", 50*scale_factor))
        .x_label_area_size(70)
        .y_label_area_size(100)
        .build_cartesian_2d((prm.w_range.0/prm.w_0)..(prm.w_range.1/prm.w_0), (min_y .. max_y).log_scale())?;

    chart.configure_mesh()
        .x_label_style(("helvetica", 50))
        .y_label_style(("helvetica", 50))
        .disable_mesh()
        .set_all_tick_mark_size(20)
        // .bold_line_style(original_style)
        .x_label_formatter(&|v| format!("{0:.2}", v))
        .y_label_formatter(&|v| format!("{:e}", v))
        .draw()?;

    let col = vec![LIGHTBLUE_A700,LIME_A700,DEEPORANGE_A400,BLACK];

    ardofs.iter().enumerate().for_each(|(ind, (qual, ardof))| {   
        let x: Vec<f64> = (omegas.to_owned() / prm.w_0).to_vec();
        let y: Vec<f64> = ardof.to_vec();
        let data: Vec<(f64,f64)> = x.into_iter().zip(y).collect();

        let sty = col[ind].stroke_width(4);

        chart.draw_series(LineSeries::new((data.clone()).into_iter(), sty)).unwrap()
            .label(format!("Quality Factor: {}", qual))
            .legend( move |(x1, y1)| PathElement::new(vec![(x1 - 100 , y1), (x1 + 20, y1)], sty));
    });
    
    chart.draw_series(DashedLineSeries::new(vec![(1.0,min_y) , (1.0, max_y)].into_iter(), 6, 10, ShapeStyle {color: BLACK.mix(1.0), filled: true, stroke_width: 2}))?;

    chart.configure_series_labels()
        .position(SeriesLabelPosition::UpperRight)
        .label_font(("helvetica", 50))
        .border_style(&WHITE.mix(0.0))
        .background_style(&WHITE.mix(0.8))
        .draw()?;
    
    Ok(())
}

fn plot_weights_quality(weights_quality:&Vec<Array1<f64>>, prm:&Parameters, omegas:&Array1<f64>, qualities:&Vec<f64>) ->  Result<(), Box<dyn std::error::Error>> {
    
    // Generate x-axis grid
    let omegas = omegas.slice(s![..;prm.n_w / prm.n_w_bins]);

    // Find maximum y-value
    let max_q_ind = qualities.iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(index, _)| index);
    // let test = weights_quality[max_q_ind.unwrap()].to_vec();
    let max_y = *weights_quality[max_q_ind.unwrap()]
        .to_vec()
        .iter()
        .max_by(|a,b| a.total_cmp(b))
        .unwrap();

    let min_y = *weights_quality[max_q_ind.unwrap()]
        .to_vec()
        .iter()
        .min_by(|a,b| a.total_cmp(b))
        .unwrap();

    let root = SVGBackend::new("weights.svg", (1440,1080)).into_drawing_area();
        
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(40)
        // .caption("Test Dispersion", ("helvetica", 50*scale_factor))
        .x_label_area_size(70)
        .y_label_area_size(100)
        .build_cartesian_2d((prm.w_range.0/prm.w_0)..(prm.w_range.1/prm.w_0), (min_y .. max_y).log_scale())?;

    chart.configure_mesh()
        .x_label_style(("helvetica", 50))
        .y_label_style(("helvetica", 50))
        .disable_mesh()
        .set_all_tick_mark_size(20)
        // .bold_line_style(original_style)
        .x_label_formatter(&|v| format!("{0:.2}", v))
        .y_label_formatter(&|v| format!("{:e}", v))
        .draw()?;

    let col = vec![LIGHTBLUE_A700,LIME_A700,DEEPORANGE_A400,BLACK];

    qualities.iter().enumerate().for_each(|(ind, qual)| {   
        let x: Vec<f64> = (omegas.to_owned() / prm.w_0).to_vec();
        let y: Vec<f64> = weights_quality[ind].to_vec();
        let data: Vec<(f64,f64)> = x.into_iter().zip(y).collect();

        let sty = col[ind].stroke_width(4);

        chart.draw_series(LineSeries::new((data.clone()).into_iter(), sty)).unwrap()
            .label(format!("Quality Factor: {}", qual))
            .legend( move |(x1, y1)| PathElement::new(vec![(x1 - 100 , y1), (x1 + 20, y1)], sty));
    });
    
    chart.draw_series(DashedLineSeries::new(vec![(1.0,min_y) , (1.0, max_y)].into_iter(), 6, 10, ShapeStyle {color: BLACK.mix(1.0), filled: true, stroke_width: 2}))?;

    chart.configure_series_labels()
        .position(SeriesLabelPosition::UpperRight)
        .label_font(("helvetica", 50))
        .border_style(&WHITE.mix(0.0))
        .background_style(&WHITE.mix(0.8))
        .draw()?;
    
    Ok(())
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
        .y_label_area_size(150)
        .build_cartesian_2d(prm.w_range.0..prm.w_range.1, 0.0 .. weights.max())?;

    chart.configure_mesh()
    .x_label_style(("helvetica", 50))
    .y_label_style(("helvetica", 50))
    .disable_mesh()
    .set_all_tick_mark_size(20)
    // .bold_line_style(original_style)
    .x_label_formatter(&|v| format!("{0:.2}", v))
    .y_label_formatter(&|v| format!("{:e}", v))
    .draw()?;

    chart.draw_series(LineSeries::new((data.clone()).into_iter(), BLUE.stroke_width(4)))?;

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
    
    let _test = weights.to_vec();
    let _test1 = omegas.to_vec();
    let _test2 = q_pars.to_vec();
    // let _test3 = weights.to_vec();

    weights
}

fn plot_disp (dispersion:Array2<f64>, omegas:Array1<f64>, q_pars: Array1<f64>) {
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
    pub w_range: (f64,f64),
    pub routine: String,
}