#[derive(Clone)]
pub struct Parameters {
    pub l_c: f64,
    pub c: f64,
    pub r: f64,
    pub w_c: f64,
    pub n_w: usize,
    pub n_q: usize,
    pub n_w_bins: usize,
    pub del_k: f64,
    pub quality: f64,
    pub q_range: (f64,f64),
    pub w_range: (f64,f64),
    pub routine: String,
    // pub coupling: Option<f64>,
    // pub w_0: Option<f64>
}

#[derive(Clone)]
pub struct VscParameters {
    pub coupling: f64,
    pub w_0: f64,
    pub k_0: f64,
    pub beta: f64,
    pub unit: Units,
}

#[derive(Clone)]
pub enum Units {
    AU,
    CM
}

impl VscParameters {
    pub fn to_au(&mut self) {
        match &self.unit {
            Units::AU => return,
            Units::CM => {
                self.beta     = cm_to_au(&self.beta);
                self.coupling = cm_to_au(&self.coupling);
                self.w_0      = cm_to_au(&self.w_0);
                self.unit     = Units::AU;
            }
        }
    }

    pub fn to_cm(&mut self) {
        match &self.unit {
            Units::CM => return,
            Units::AU => {
                self.beta     = au_to_cm(&self.beta);
                self.coupling = au_to_cm(&self.coupling);
                self.w_0      = au_to_cm(&self.w_0);
                self.unit     = Units::CM;
            }
        }
    }
}

fn cm_to_au(val:&f64) -> f64 {
    val / 219474.63
}

fn au_to_cm(val:&f64) -> f64 {
    val * 219474.63
}