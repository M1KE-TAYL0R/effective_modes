mod backend;

mod ell_n;
use ell_n::{single_q_factor,many_q_factors};

mod dispersion_plot;
use dispersion_plot::plot_disp_wrapper;

mod ardof;
use ardof::ardof;

mod vsc_rate;
use vsc_rate::vsc_rate;

mod parameters;
use parameters::{Parameters, VscParameters,Units};

fn define_parameters() -> Parameters {
    let prm = Parameters{
        l_c: 0.0,                         // Placeholder
        c: 1.0 / 137.0,                   // Speed of light in atomic units
        r: 0.0,                           // Placeholder if quality != 0
        // w_c: 0.1,                         // Cavity resonance frequency in atomic units
        w_c: 1100. / 219474.63,
        // n_w: 2000,                       // Number of \omega grid points
        // n_q: 2000,                      // Number of q_\| grid points per bin integration
        // n_w_bins: 2000,                   // Number of omega_n bins
        n_w: 2000,                       // Number of \omega grid points
        n_q: 240,                      // Number of q_\| grid points per bin integration
        n_w_bins: 2000,                   // Number of omega_n bins
        del_k: 0.0,                       // Placeholder
        quality: 50.0,                   // Cavity Quality Factor
        q_range: (0.0,0.5),              // Range of q_\| points integrated over
        w_range: (100. / 219474.63, 500. / 219474.63),            // Range of omega_n
        // q_range: (-5.0,5.0),              // Range of q_\| points integrated over
        // w_range: (0.099, 0.108),            // Range of omega_n
        // q_range: (0.0,10.0),              // Range of q_\| points integrated over
        // w_range: (0.09, 0.12),            // Range of omega_n
        routine: "VSC_k".to_string(),      // Set the routine: ManyQ, SingQ, Dispn, ARDOF, VSC_k
        // coupling: None,
        // w_0: None
    };

    prm
}

fn define_vsc_parameters() -> VscParameters {
    let vsc_prm = VscParameters{
        coupling: 0.0,
        w_0: 1200.,
        k_0: 5.946954192406803128e-08, // From Wenxiang HEOM
        beta: 200.0,
        unit: Units::CM,
        // k_vsc: Vec::new()
    };

    vsc_prm
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
        "VSC_k"  => vsc_rate(prm,define_vsc_parameters()),
        _        => panic!("Invalid Routine Parameter: check prm.routine")
    };
}

