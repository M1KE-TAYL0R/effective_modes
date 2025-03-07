use ndarray::Array1;
use crate::{parameters::Parameters, backend::enhancement_function};
use std::f64::consts::PI;
use rayon::prelude::*;
use plotters::{prelude::*, style::{full_palette::{DEEPORANGE_A400, LIGHTBLUE_A700, LIME_A700}, Color}};

pub fn ardof(mut prm: Parameters) {
    let q_0 = prm.w_0/prm.c;
    prm.l_c  = 2.0 * PI / q_0;

    let omegas = Array1::linspace(prm.w_range.0, prm.w_range.1,prm.n_w);
    let q_pars = Array1::linspace(prm.q_range.0, prm.q_range.1, prm.n_q);

    let qualities = vec![5000.,200.,100.,50.];

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
                
                let jacobian =  q_par * omega / prm.c / q_z * dq;

                *temp1 = enhancement_function(*omega, q_par, &prm) * jacobian;
            });
            ardofs[ind].1[ind1] =  temp.iter().sum();
            // ardofs[ind].1[ind1] *= dq;

        });

        // ardofs[ind].1 *= dq;
    });

    plot_ardofs_quality(&ardofs, &prm, &omegas).unwrap();
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

    let scale: u32 = 10;
    
    let root = SVGBackend::new("ardofs.svg", (scale*1440,scale*1080)).into_drawing_area();

    let original_style = ShapeStyle {
        color: BLACK.into(),
        filled: false,
        stroke_width: 2*scale,
    };
    
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(40*scale)
        // .caption("Test Dispersion", ("helvetica", 50*scale_factor))
        .x_label_area_size(70*scale)
        .y_label_area_size(100*scale)
        .build_cartesian_2d((prm.w_range.0/prm.w_0)..(prm.w_range.1/prm.w_0), (min_y .. max_y).log_scale())?;

    chart.configure_mesh()
        .x_label_style(("helvetica", 50*scale))
        .y_label_style(("helvetica", 50*scale))
        .disable_mesh()
        .set_all_tick_mark_size(20*scale)
        // .bold_line_style(original_style)
        // .light_line_style(original_style)
        .axis_style(original_style)
        .x_label_formatter(&|v| format!("{0:.2}", v))
        .y_label_formatter(&|v| format!("{:e}", v))
        .draw()?;

    let col = vec![LIGHTBLUE_A700,LIME_A700,DEEPORANGE_A400,BLACK];

    ardofs.iter().enumerate().for_each(|(ind, (qual, ardof))| {   
        let x: Vec<f64> = (omegas.to_owned() / prm.w_0).to_vec();
        let y: Vec<f64> = ardof.to_vec();
        let data: Vec<(f64,f64)> = x.into_iter().zip(y).collect();

        let sty = col[ind].stroke_width(4*scale);

        chart.draw_series(LineSeries::new((data.clone()).into_iter(), sty)).unwrap()
            .label(format!("Quality Factor: {}", qual))
            .legend( move |(x1, y1)| PathElement::new(vec![(x1 - 100*(scale as i32) , y1), (x1 - 20 * scale as i32, y1)], sty));
    });
    
    chart.draw_series(DashedLineSeries::new(vec![(1.0,min_y) , (1.0, max_y)].into_iter(), 
        6*scale, 
        10*scale, 
        ShapeStyle {color: BLACK.mix(1.0), 
        filled: true, stroke_width: 2*scale}))?;

    chart.configure_series_labels()
        .position(SeriesLabelPosition::Coordinate(850 * scale as i32, 700 * scale as i32))
        .margin(10 * scale)
        .label_font(("helvetica", 50*scale))
        .border_style(&WHITE.mix(0.0))
        .background_style(&WHITE.mix(0.0))
        .draw()?;
    
    Ok(())
}