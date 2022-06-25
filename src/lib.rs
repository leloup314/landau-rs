// Bring probability-density functions into the namespace
mod pdfs;

// Make PDFs public
pub use crate::pdfs::{gauss_pdf, landau_pdf, langau_pdf};

/// Check for correct parameters
fn check_parameters(eta: &mut f64, sigma: &mut f64, amplitude: Option<f64>) {

    // Check eta
    if *eta < 1e-9 {
        println!("WARNING: eta < 1e-9 is not supported. Setting eta = 1e-9.");
        *eta = 1e-9;
    }
    // Check sigma
    if *sigma < 0.0 {
        *sigma *= -1.0;
    }
    if *sigma > 100.0 * *eta {
        println!("WARNING: sigma > 100 * eta can lead to oszillations. Check result.")
    }
    // Check amplitude
    let amp_check = amplitude.unwrap_or(1.0);
    if amp_check < 0.0 {
        panic!("Amplitude has to be >= 0, got {}", amp_check);
    }
}

fn scale_to_mpv(mu: f64, eta: f64, sigma: f64, amplitude: Option<f64>) -> f64 {
    1 as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
        let a = 15f64;
        let b = 1f64;
        let c = 1f64;
        let d = 1f64;
        assert_eq!(langau_pdf(a, b, c, d), 0.006297828887560015);
    }
}
