use crate::{parameters::Parameters, backend::{set_cavity_params,calc_dispersion,bin_dispersion,par_weights_gen}};
use std::f64::consts::PI;
use ndarray::{Array1,s};
use ndarray_npy::{write_npy,read_npy};
use std::path::Path;
use plotters::{prelude::*, style::{full_palette::{DEEPORANGE_A400, LIGHTBLUE_A700, LIME_A700}, Color}};
use statrs::statistics::Statistics;


pub fn many_q_factors(mut prm:Parameters) {

    prm.w_range = (prm.w_c * 0.9, prm.w_c * 1.4);

    let q_0 = prm.w_c/prm.c;
    prm.l_c  = 2.0 * PI / q_0;

    let omegas = Array1::linspace(prm.w_range.0, prm.w_range.1,prm.n_w);

    let qualities = vec![5000.,500.,50.];

    let (prm,output) = scan_quality_factors(&qualities, prm, &omegas);

    save_output(&prm,output);
}

fn save_output (prm: &Parameters, output: Vec<(f64,Array1<f64>)>) {

    output.iter().for_each(|(qual,weight)| {
        let fname = format!("data/ell_n__w{}_nq{}_qual{}.npy", prm.w_c, prm.n_q, qual);
        println!("Saving: '{}'", fname);

        write_npy(fname, weight).unwrap();
    });
}

pub fn single_q_factor(mut prm:Parameters) {

    prm = set_cavity_params(prm);

    let omegas = Array1::linspace(prm.w_range.0, prm.w_range.1,prm.n_w);
    let q_pars = Array1::linspace(prm.q_range.0, prm.q_range.1, prm.n_q);

    let dispersion = calc_dispersion(&prm, &omegas, &q_pars);

    let weights = bin_dispersion(&prm, &omegas, &q_pars, &dispersion);
    plot_weights(weights, &prm, &omegas).unwrap();

}

pub fn scan_quality_factors(qualities:&Vec<f64>, mut prm:Parameters, omegas:&Array1<f64>) -> (Parameters, Vec<(f64,Array1<f64>)>){

    let mut weights_quality:Vec<Array1<f64>> = Vec::new();

    qualities.iter().for_each(|quality| {
        prm.quality = *quality;

        prm = set_cavity_params(prm.clone());

        // weights_quality.push(weights_gen(&prm, omegas, q_pars));

        weights_quality.push(par_weights_gen(&prm, omegas));

        // println!("Calculating Q = {}", quality);

        // let fname = format!("ell_n__w{}_nq{}_qual{}.npy", prm.w_c, prm.n_q, prm.quality);

        // if Path::new(fname.as_str()).is_file() {
        //     let weights: Array1<f64> = read_npy(fname).unwrap();
        //     weights_quality.push(weights);
        // }
        // else {
        //      // Define r
        //      let gamma = prm.w_c / prm.quality;
        //      prm.del_k = gamma / prm.c;
 
        //      // Lorentzian linewidth approx
        //      prm.r = (- 2.0 * PI / prm.quality).exp().powf(0.25);

        //     let dispersion = calc_dispersion(&prm, &omegas, &q_pars);
 
        //     weights_quality.push(bin_dispersion(&prm, &omegas, &q_pars, &dispersion));
        // }
    });

    plot_weights_quality(&weights_quality, &prm, &omegas, &qualities).unwrap();

    (prm, qualities.clone().into_iter().zip(weights_quality).collect())

}

pub fn _weights_gen(prm: &Parameters, omegas:&Array1<f64>, q_pars: &Array1<f64>) -> Array1<f64>{
    println!("Calculating Q = {}", prm.quality);

    let fname = format!("ell_n__w{}_nq{}_qual{}.npy", prm.w_c, prm.n_q, prm.quality);

    if Path::new(fname.as_str()).is_file() {
        let weights: Array1<f64> = read_npy(fname).unwrap();
        return weights
    }
    else {

        let dispersion = calc_dispersion(&prm, omegas, q_pars);

        return bin_dispersion(&prm, omegas,&q_pars, &dispersion);
    }
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

    println!("Max ell_n = {}", {max_y});

    let min_y = *weights_quality[max_q_ind.unwrap()]
        .to_vec()
        .iter()
        .min_by(|a,b| a.total_cmp(b))
        .unwrap();

    let scale: u32 = 10;

    let original_style = ShapeStyle {
        color: BLACK.into(),
        filled: false,
        stroke_width: 4*scale,
    };

    let root = SVGBackend::new("weights.svg", (1440*scale,1080*scale)).into_drawing_area();
        
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(70*scale)
        // .caption("Test Dispersion", ("helvetica", 50*scale_factor))
        .x_label_area_size(70*scale)
        .y_label_area_size(100*scale)
        .build_cartesian_2d((prm.w_range.0/prm.w_c)..(prm.w_range.1/prm.w_c), (min_y .. max_y).log_scale())?;

    chart.configure_mesh()
        .x_label_style(("helvetica", 70*scale))
        .y_label_style(("helvetica", 70*scale))
        .disable_mesh()
        .set_all_tick_mark_size(30*scale)
        .axis_style(original_style)
        .x_label_formatter(&|v| format!("{0:.2}", v))
        .y_label_formatter(&|v| format!("{:e}", v))
        .draw()?;

    let col = vec![LIGHTBLUE_A700,LIME_A700,DEEPORANGE_A400,BLACK];

    qualities.iter().enumerate().for_each(|(ind, qual)| {   
        let x: Vec<f64> = (omegas.to_owned() / prm.w_c).to_vec();
        let y: Vec<f64> = weights_quality[ind].to_vec();
        let data: Vec<(f64,f64)> = x.into_iter().zip(y).collect();

        let sty = col[ind].stroke_width(8*scale);

        chart.draw_series(LineSeries::new((data.clone()).into_iter(), sty)).unwrap()
            .label(format!("Quality Factor: {}", qual))
            .legend( move |(x1, y1)| PathElement::new(vec![(x1 - 100*scale as i32 , y1), (x1 - 20*scale as i32, y1)], sty));
    });
    
    chart.draw_series(DashedLineSeries::new(vec![(1.0,min_y) , (1.0, max_y)].into_iter(), 6*scale, 10*scale, ShapeStyle {color: BLACK.mix(1.0), filled: true, stroke_width: 2*scale}))?;

    chart.configure_series_labels()
        .position(SeriesLabelPosition::UpperRight)
        .label_font(("helvetica", 70*scale))
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
