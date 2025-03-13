use crate::backend::set_cavity_params;
use crate::ell_n::scan_quality_factors;
use crate::parameters::Parameters;

use std::f64::consts::PI;
use ndarray::Array1;

pub fn vsc_rate(mut prm: Parameters) {
    let q_0 = prm.w_c/prm.c;
    prm.l_c  = 2.0 * PI / q_0;

    let omegas = Array1::linspace(prm.w_range.0, prm.w_range.1,prm.n_w);
    let q_pars = Array1::linspace(prm.q_range.0, prm.q_range.1, prm.n_q);

    let qualities = vec![5000.,500.,50.];
    let couplings = vec![3.125e-4, 6.25e-4, 9.375e-4, 1.25e-3];

    let (mut prm,ell_n_quality) = scan_quality_factors(qualities, prm, &omegas, q_pars);

    let mut k_vsc_quality: Vec<(f64,f64,Array1<f64>)> = Vec::new(); // (coupling, quality, ell_n)

    let w_0 = prm.w_0.unwrap();

    let shape = (*(&omegas).shape())[0];

    couplings.iter().for_each(|coupling| {
        prm.coupling = Some(*coupling);
        ell_n_quality.iter().for_each(|(quality, ell_n)| {

            // let t_c = - prm.c / 2.0 / prm.l_c / (prm.r.powi(2)).ln();
            
            let k_vsc: Array1<f64> = Array1::zeros(shape);

            omegas.clone().iter().for_each(|w_0| {
                
            });

            k_vsc_quality.push((*coupling,*quality, ell_n.clone()));
        });
    });
    

}

fn _calc_k_vsc(omegas:&Array1<f64>, w_0: f64, ell_n: &Array1<f64>, quality:f64, mut prm:Parameters) {

    prm.quality = quality;
    prm = set_cavity_params(prm);


    
}

fn calc_t_c(prm:&Parameters) -> f64 {
    
    let t_c = - prm.c / 2.0 / prm.l_c / (prm.r.powi(2)).ln();

    t_c
}