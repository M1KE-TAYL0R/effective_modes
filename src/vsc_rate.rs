use crate::backend::{set_cavity_params,par_weights_gen};
// use crate::ell_n::weights_gen;
use crate::parameters::{Parameters,VscParameters};
use plotters::{prelude::*, style::{full_palette::{DEEPORANGE_A400, LIGHTBLUE_A700, LIME_A700}, Color}};
use indicatif::{ProgressIterator, ProgressStyle};

use std::collections::HashMap;
use std::f64::consts::PI;
use ndarray::Array1;


pub fn vsc_rate(mut prm: Parameters, mut vsc_prm: VscParameters) {
    let q_0 = prm.w_c/prm.c;
    prm.l_c  = 2.0 * PI / q_0;
    
    // let qualities = vec![500.,200.,50.];
    // let couplings = vec![3.125e-4, 6.25e-4, 9.375e-4, 1.25e-3];

    let qualities = vec![500.,200., 50.];
    let couplings = vec![1.25e-3];

    vsc_prm.to_au();

    let omegas = Array1::linspace(prm.w_range.0, prm.w_range.1,prm.n_w_bins);
    let q_pars = Array1::linspace(prm.q_range.0, prm.q_range.1, prm.n_q);

    // let mut k_vsc_qual_coup: Vec<(f64,f64,Array1<f64>)> = Vec::with_capacity(qualities.len() * couplings.len()); // (coupling, quality, k_vsc(w_c) )
    // let mut k_vsc_qual_coup: Vec<(f64,f64,Array1<f64>)> = Vec::new();
    let mut k_vsc_qual_coup: HashMap<(usize,usize), Array1<f64>> = HashMap::new();

    // Iterate over couplings, w_c, quality
    qualities.iter().enumerate().for_each(|(qual_ind, quality)| {
        prm.quality = *quality;

        couplings.iter().enumerate().for_each(|(coup_ind, coupling)| {
            vsc_prm.coupling = *coupling;

            let mut k: Array1<f64> = Array1::zeros(omegas.len());

            
            let prog_style = ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})").unwrap();

            omegas.clone().indexed_iter().progress_with_style(prog_style).for_each(|(w_ind,w_c)| {
                prm.w_c = *w_c;
                prm.w_range = (prm.w_c * 0.9, prm.w_c * 1.4);
                
                prm = set_cavity_params(prm.clone());

                let w_scan = Array1::linspace(prm.w_range.0, prm.w_range.1,prm.n_w);          

                let ell_n = par_weights_gen(&prm, &w_scan);

                let _temp = q_pars.to_vec();
                let _temp2 = ell_n.to_vec();
                let _temp1: f64 =  ell_n.to_vec().iter().sum();

                k[w_ind] = calc_k(&w_scan,&vsc_prm,&ell_n,&prm);

                if k[w_ind].is_nan() {
                    println!("l_n/ V = {}", ell_n);
                    println!("omega = {}", w_c);
                }
            });
            
            k_vsc_qual_coup.insert((coup_ind,qual_ind),k);
        });
    });

    vsc_prm.to_cm();
    let omegas_cm = omegas * 219474.63;

    plot_k_vsc(&k_vsc_qual_coup,&omegas_cm,qualities,couplings).unwrap();

}

fn plot_k_vsc(k_vsc_qual_coup: &HashMap<(usize,usize), Array1<f64>> , omegas_cm: &Array1<f64>, qualities: Vec<f64>, _couplings: Vec<f64>) -> Result<(), Box<dyn std::error::Error>>{

    // Plot k_vsc vs. quality factor
    let coup_ind = 0;

    // Filter out k_vsc
    let mut k_vsc_qual = k_vsc_qual_coup.clone();
    k_vsc_qual.retain(|key, _| key.0 == coup_ind );

    // Find maximum y-value
    let maxes_y: Vec<((usize,usize),f64)> = k_vsc_qual
        .iter()
        .map(|(a,b)| 
            (*a,*b.to_vec().iter().max_by(|a,b| a.total_cmp(b))
            .unwrap()))
            .collect();

    let max_y = *maxes_y
        .iter()
        .max_by(|(_,a),(_,b)| a.total_cmp(b))
        .unwrap();

    println!("Max k_vsc = {}", {max_y.1});

    let scale: u32 = 10;

    let original_style = ShapeStyle {
        color: BLACK.into(),
        filled: false,
        stroke_width: 4*scale,
    };

    let root = SVGBackend::new("k_qual.svg", (1440*scale,1080*scale)).into_drawing_area();
    
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(70*scale)
        // .caption("Test Dispersion", ("helvetica", 50*scale_factor))
        .x_label_area_size(70*scale)
        .y_label_area_size(100*scale)
        .build_cartesian_2d((omegas_cm[0])..*(omegas_cm.last().unwrap()), (0.01 .. max_y.1).log_scale())?;

    chart.configure_mesh()
        .x_label_style(("helvetica", 20*scale))
        .y_label_style(("helvetica", 20*scale))
        .disable_mesh()
        .set_all_tick_mark_size(30*scale)
        .axis_style(original_style)
        .x_label_formatter(&|v| format!("{0:.2}", v))
        .y_label_formatter(&|v| format!("{:e}", v))
        .draw()?;

    let col = vec![LIGHTBLUE_A700,LIME_A700,DEEPORANGE_A400,BLACK];

    
    qualities.iter().enumerate().for_each(|(ind, qual)| {
        let x: Vec<f64> = omegas_cm.to_owned().to_vec();
        let y: Vec<f64> = k_vsc_qual_coup[&(coup_ind,ind)].to_vec();
        let data: Vec<(f64,f64)> = x.into_iter().zip(y).collect();

        let sty = col[ind].stroke_width(8*scale);

        chart.draw_series(LineSeries::new((data.clone()).into_iter(), sty)).unwrap()
            .label(format!("Quality Factor: {}", qual))
            .legend( move |(x1, y1)| PathElement::new(vec![(x1 - 100*scale as i32 , y1), (x1 - 20*scale as i32, y1)], sty));
    });

    Ok(())
}


fn calc_k(omegas:&Array1<f64>, vsc_prm: &VscParameters, ell_n: &Array1<f64>, prm:&Parameters) -> f64 {

    let mut k = 0.0;

    let t_c = calc_t_c(&prm);
    let w_0 = &vsc_prm.w_0;
    let beta = &vsc_prm.beta;

    let dw = (prm.w_range.1 - prm.w_range.0) / omegas.len() as f64;

    omegas.indexed_iter().for_each(|(ind,omega)|{
        let _t1 = lorentzian(omega, w_0, t_c);
        let _t2 = bose_einstein_dist(beta, w_0);
        // let temp = t1 * t2 ;
        // if temp.is_nan(){
        //     println!("temp is NaN");
        // }
        // k += ell_n[ind] * lorentzian(omega, w_0, t_c) * bose_einstein_dist(beta, w_0);

        let _temp = ell_n.to_vec();
        let _temp2 = ell_n[ind];

        k += ell_n[ind] * lorentzian(omega, w_0, t_c) * bose_einstein_dist(beta, w_0);// * (- beta * omega).exp();
    });

    // println!("k = {}",k);

    // k * 4.0 * vsc_prm.coupling.powi(2) * dw / vsc_prm.k_0
    1.0 + k * 4.0 * vsc_prm.coupling.powi(2)* dw / vsc_prm.k_0
}

fn calc_t_c(prm:&Parameters) -> f64 {
    
    let t_c = -  prm.l_c / prm.c / (prm.r.powi(2)).ln();

    t_c
}

fn lorentzian(omega:&f64, w_0: &f64, t_c:f64) -> f64 {
    let numerator = omega*omega * w_0 / t_c;
    let denominator = (omega*omega - w_0*w_0).powi(2) + (w_0/t_c).powi(2);

    numerator/denominator
}

fn bose_einstein_dist(beta:&f64, w_0: &f64) -> f64{
    let n = ((beta * w_0).exp() - 1.0).powi(-1);

    n
}