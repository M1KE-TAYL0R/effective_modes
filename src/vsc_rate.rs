use crate::backend::{set_cavity_params,par_weights_gen};
// use crate::ell_n::weights_gen;
use crate::parameters::{Parameters,VscParameters};
use plotters::{prelude::*, style::{full_palette::{DEEPORANGE_A400, LIGHTBLUE_A700, LIME_A700}, Color}};

use std::f64::consts::PI;
use ndarray::Array1;

pub fn vsc_rate(mut prm: Parameters, mut vsc_prm: VscParameters) {
    let q_0 = prm.w_c/prm.c;
    prm.l_c  = 2.0 * PI / q_0;
    
    let qualities = vec![5000.,500.,50.];
    let couplings = vec![3.125e-4, 6.25e-4, 9.375e-4, 1.25e-3];

    vsc_prm.to_au();

    let omegas = Array1::linspace(prm.w_range.0, prm.w_range.1,prm.n_w);
    let q_pars = Array1::linspace(prm.q_range.0, prm.q_range.1, prm.n_q);

    // let mut k_vsc_qual_coup: Vec<(f64,f64,Array1<f64>)> = Vec::with_capacity(qualities.len() * couplings.len()); // (coupling, quality, k_vsc(w_c) )
    let mut k_vsc_qual_coup: Vec<(f64,f64,Array1<f64>)> = Vec::new();

    // Iterate over couplings, w_c, quality
    qualities.iter().for_each(|quality| {
        prm.quality = *quality;

        couplings.iter().for_each(|coupling| {
            vsc_prm.coupling = *coupling;

            let mut k: Array1<f64> = Array1::zeros(omegas.len());

            omegas.clone().indexed_iter().for_each(|(w_ind,w_c)| {
                prm.w_c = *w_c;          
                prm = set_cavity_params(prm.clone());

                let ell_n = par_weights_gen(&prm, &omegas, &q_pars);

                k[w_ind] = calc_k(&omegas,&vsc_prm,&ell_n,&prm);
            });
            
            k_vsc_qual_coup.push((*coupling,*quality,k));
        });
    });

    vsc_prm.to_cm();
    let omegas_cm = omegas * 219474.63;

    plot_k_vsc(&k_vsc_qual_coup,omegas_cm,qualities,couplings).unwrap();

}

fn plot_k_vsc(k_vsc_qual_coup: &Vec<(f64,f64,Array1<f64>)>, omegas_cm: Array1<f64>, qualities: Vec<f64>, couplings: Vec<f64>) -> Result<(), Box<dyn std::error::Error>>{

    // Plot k_vsc vs. quality factor
    let coup_ind = 2;

    // Filter out k_vsc
    let k_vsc_qual = k_vsc_qual_coup.iter().filter_map(|(coup,qual,k)| {
        if coup == &couplings[coup_ind] {
            Some((qual,k))
        }
        else {
            None
        }
        }).collect();

    // Find maximum y-value
    let maxes_y: Vec<(f64,f64,f64)> = *k_vsc_qual_coup
        .iter()
        .map(|(a,b,c)| 
            (*a,*b,*c.to_vec().iter().max_by(|a,b| a.total_cmp(b))
            .unwrap()))
            .collect();

    let max_y = *maxes_y
        .iter()
        .max_by(|(_,_,a),(_,_,b)| a.total_cmp(b))
        .unwrap();

    println!("Max k_vsc = {}", {max_y.2});

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
        .build_cartesian_2d((omegas_cm[0])..*(omegas_cm.last().unwrap()), 0.9 .. max_y.2)?;

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

    qualities.iter().enumerate().for_each(|(q_ind, _)| {

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
        k += ell_n[ind] * lorentzian(omega, w_0, t_c) * bose_einstein_dist(beta, w_0);
    });

    k * 4.0 * vsc_prm.coupling.powi(2) * dw / vsc_prm.k_0
}

fn calc_t_c(prm:&Parameters) -> f64 {
    
    let t_c = - prm.c / 2.0 / prm.l_c / (prm.r.powi(2)).ln();

    t_c
}

fn lorentzian(omega:&f64, w_0: &f64, t_c:f64) -> f64 {
    let numerator = omega*omega * w_0 / t_c;
    let denominator = (omega*omega - w_0*w_0).powi(2) + (w_0/t_c).powi(2);

    numerator/denominator
}

fn bose_einstein_dist(beta:&f64, w_0: &f64) -> f64{
    let n = ((beta * w_0).exp() -1.0).powi(-1);

    n
}