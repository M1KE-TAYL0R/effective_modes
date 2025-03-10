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
use parameters::Parameters;

fn define_parameters() -> Parameters {
    let prm = Parameters{
        l_c: 0.0,                         // Placeholder
        c: 1.0 / 137.0,                   // Speed of light in atomic units
        r: 0.0,                           // Placeholder if quality != 0
        w_c: 0.1,                         // Cavity resonance frequency in atomic units
        n_w: 10000,                       // Number of \omega grid points
        n_q: 100000,                      // Number of q_\| grid points per bin integration
        n_w_bins: 5000,                   // Number of omega_n bins
        del_k: 0.0,                       // Placeholder
        quality: 50.0,                   // Cavity Quality Factor
        // q_range: (-5.0,5.0),              // Range of q_\| points integrated over
        // w_range: (0.099, 0.108),            // Range of omega_n
        q_range: (0.0,10.0),              // Range of q_\| points integrated over
        w_range: (0.09, 0.12),            // Range of omega_n
        routine: "ARDOF".to_string(),      // Set the routine: ManyQ, SingQ, Dispn, ARDOF
        coupling: None,
        w_0: None
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
        "VSC_k"  => vsc_rate(prm),
        _        => panic!("Invalid Routine Parameter: check prm.routine")
    };
}

