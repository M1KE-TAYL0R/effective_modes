use std::f64::consts::PI;
use crate::parameters::Parameters;
use ndarray_linalg::c64;
use ndarray::{s, Array1, Array2};
use ndarray_npy::{write_npy,read_npy};
use std::path::Path;
use rayon::prelude::*;

pub fn set_cavity_params(mut prm:Parameters) -> Parameters {
    // Calculate prm.del_k and prm.r from prm.quality
    if prm.quality == 0.0 {
        // Lorentzian linewidth approx
        let gamma = - prm.w_c / (2.0 * PI) * prm.r.powi(4).ln(); 

        prm.del_k = gamma / prm.c;
        prm.quality = prm.w_c / gamma;
    }
    // Calculate prm.del_k and prm.quality from prm.r 
    else {
        let gamma = prm.w_c / prm.quality;
        prm.del_k = gamma / prm.c;

        // Lorentzian linewidth approx
        prm.r = (- 2.0 * PI / prm.quality).exp().powf(0.25);
    }

    // println!("Reflectance: {}", prm.r);
    // println!("Quality Factor: {}", prm.quality);
    // println!("Delta q_perp: {}", prm.del_k);

    let q_0 = prm.w_c/prm.c;
    prm.l_c  = 2.0 * PI / q_0;

    prm
}

pub fn enhancement_function (omega:f64, q_par:f64, prm:&Parameters) -> f64 {

    let q_z = (omega.powi(2) / prm.c.powi(2) - q_par.powi(2)).sqrt();

    if q_z <  PI / prm.l_c || q_z.is_nan(){
        return 0.0;
    }

    let r = prm.r;

    let mut weight: c64 = 1.0 + r * (prm.l_c * q_z * c64::i()).exp(); 
    weight /= 1.0 - r*r*(2.0 * prm.l_c * q_z * c64::i()).exp();

    let re_weight = weight.norm() - 1.0;
    // let re_weight = weight.norm();

    if re_weight < 0.0 { return 0.0}

    re_weight * re_weight
}

fn integrate_bin (mut bin:Array2<f64>, omegas_bin:Array1<f64>, q_pars: &Array1<f64>, prm:&Parameters) -> f64 {
    
    let dq = (prm.q_range.1 - prm.q_range.0) / (prm.n_q as f64);
    let dw = (prm.w_range.1 - prm.w_range.0) / (prm.n_w as f64);

    bin.indexed_iter_mut().for_each(|(ind, val)| {
        let q_par = q_pars[ind.1];
        let omega = omegas_bin[ind.0];
        let q_z = (omega.powi(2) - prm.c*prm.c * q_par.powi(2)).sqrt();
        
        let mut jacobian =  q_par * omega / prm.c / q_z; //.powf(-1.0);

        
        let qn_2 = (omega / prm.c).powf(2.0);
        let qn_z_2 = (2.0 * PI / prm.l_c).powf(2.0);
        let qn_par = (qn_2 - qn_z_2).sqrt();
        
        if prm.del_k < qn_par {

            let del_theta = 2.0 * (prm.del_k / qn_par).asin();

            jacobian *= del_theta;
        }
        else {
            jacobian *= PI;
        }

        *val *= jacobian;
    });

   let weight = bin.sum() * dq * dw;

   weight
}

pub fn bin_dispersion (prm:&Parameters, omegas:&Array1<f64>, q_pars: &Array1<f64>, dispersion:&Array2<f64>) -> Array1<f64>{
    
    let bin_size = prm.n_w / prm.n_w_bins;
    let mut weights: Array1<f64> = Array1::zeros( prm.n_w_bins);

    for bin_ind in 0..prm.n_w_bins {
        let omegas_bin = omegas.slice(s![(bin_ind * bin_size) .. ((bin_ind+1)*bin_size)]);
        let bin = dispersion.slice(s![(bin_ind * bin_size) .. ((bin_ind+1)*bin_size),..]);

        weights[bin_ind] = integrate_bin(bin.to_owned(), omegas_bin.to_owned(), &q_pars, prm);
    }
    
    weights
}

pub fn calc_dispersion (prm:&Parameters, omegas:&Array1<f64>, q_pars: &Array1<f64>) -> Array2<f64>{
    let mut disp: Array2<f64> = Array2::zeros((prm.n_q,prm.n_w));

    disp.indexed_iter_mut().for_each(|(ind, val)| {
        let q_par = q_pars[ind.0];
        let omega = omegas[ind.1];

        *val = enhancement_function(omega, q_par, &prm);
    });

    disp.t().to_owned()
}

pub fn par_weights_gen(prm: &Parameters, omegas:&Array1<f64>, q_pars: &Array1<f64>) -> Array1<f64>{
    // println!("Calculating Q = {}", prm.quality);
    let mut q_pars = q_pars.clone();
    let num_std = 5.0;
    let read_files = true;

    let fname = format!("data2/ell_n__w{}_nq{}_qual{}.npy", prm.w_c, prm.n_q, prm.quality);

    if Path::new(fname.as_str()).is_file() && read_files{
        // println!("Reading: '{}'", fname);

        let weights: Array1<f64> = read_npy(fname).unwrap();
        return weights
    }
    else {
        // std::env::set_var("RAYON_NUM_THREADS", "96");
        let mut dq = (prm.q_range.1 - prm.q_range.0) / (prm.n_q as f64);
        let dw = (prm.w_range.1 - prm.w_range.0) / (prm.n_w as f64);

        let bin_size = prm.n_w / prm.n_w_bins;
        let bins: Array1<f64> = Array1::zeros(prm.n_w_bins);

        let weights = bins.iter().enumerate().map(|(bin_ind,_)| {

            let mut weight: f64 = 0.0;
            
            let omegas_bin = omegas.slice(s![(bin_ind * bin_size) .. ((bin_ind+1)*bin_size)]);

            omegas_bin.iter().for_each(| omega|  {

                let mut integrands: Vec<f64> = Array1::zeros(prm.n_q).to_vec();

                // Define range of q_pars:
                let w_bound = num_std * prm.c * prm.del_k;
                // let w_bound = prm.w_c * 0.2;
                if *omega > prm.w_c + w_bound {
                    let q_plus = ((omega + 3.0 *  w_bound).powi(2) - prm.w_c.powi(2)).sqrt() / prm.c;
                    let q_minus = ((omega - w_bound).powi(2) - prm.w_c.powi(2)).sqrt() / prm.c;
                    q_pars = Array1::linspace(q_minus, q_plus, prm.n_q);
                    let _temp = w_bound / prm.w_c;
                    dq = (q_plus - q_minus)/ (prm.n_q as f64);
                }
                else if *omega > prm.w_c - w_bound {
                    // // q_pars = Array1::linspace(0.0, 0.5, prm.n_q);
                    // q_pars = Array1::linspace(0.0, num_std* prm.del_k, prm.n_q);
                    // dq = (num_std* prm.del_k)/ (prm.n_q as f64);

                    let q_plus = ((omega + 3.0* w_bound).powi(2) - prm.w_c.powi(2)).sqrt() / prm.c;
                    q_pars = Array1::linspace(0.0, q_plus, prm.n_q);
                    dq = q_plus / (prm.n_q as f64);
                }
                else {
                    // let q_plus = ((prm.w_c +  w_bound).powi(2) - prm.w_c.powi(2)).sqrt() / prm.c;
                    // q_pars = Array1::linspace(0.0, q_plus, prm.n_q);
                    dq = 0.5 / (prm.n_q as f64);
                    q_pars = Array1::linspace(0.0, 0.5, prm.n_q);
                }

                q_pars.to_vec().into_par_iter().zip(integrands.par_iter_mut()).for_each(|(q_par,integrand)| {
    
                    let q_z = (omega.powi(2) - prm.c*prm.c * q_par.powi(2)).sqrt();
                    
                    let mut jacobian =  q_par * omega / prm.c / q_z;
                    if q_z.is_nan() { jacobian = 0.0}
    
                    let qn_2 = (omega / prm.c).powi(2);
                    let qn_z_2 = (2.0 * PI / prm.l_c).powi(2);
                    let qn_par = (qn_2 - qn_z_2).sqrt();
                    
                    if prm.del_k < qn_par {
    
                        let del_theta = 2.0 * (prm.del_k / qn_par).asin();
    
                        jacobian *= del_theta;
                    }
                    else {
                        jacobian *= PI;
                    }
                    let temp = enhancement_function(*omega, q_par, &prm) * jacobian;
                    *integrand = temp;
                });
                let _temp = integrands.to_vec();
                let sum = integrands.iter().sum::<f64>();
                weight +=  sum * dq;
            });

            weight
        }).collect::<Array1<f64>>();

        let _temp = weights.to_vec();

        let out = weights  * dw;

        // println!("Saving: '{}'", fname);

        write_npy(fname, &out).unwrap();

        out
    }
}