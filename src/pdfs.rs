/*
This file contains the probability-density-functions (PDFs) used within landau-rs
Taken from LCG ROOT MathLib
License info:
Authors: Andras Zsenei & Lorenzo Moneta   06/2005

**********************************************************************
*                                                                    *
* Copyright (c) 2005 , LCG ROOT MathLib Team                         *
*                                                                    *
*                                                                    *
**********************************************************************

Langaus authors:
Based on a Fortran code by R.Fruehwirth (fruhwirth@hephy.oeaw.ac.at)
Adapted for C++/ROOT by H.Pernegger (Heinz.Pernegger@cern.ch) and
Markus Friedl (Markus.Friedl@cern.ch)

Adaption for Python by David-Leon Pohl, pohl@physik.uni-bonn.de
Adaption for Rust by Pascal Wolf, wolf@physik.uni-bonn.de
*/

// Rust thinks the constants as well as pub functions in this crate are unused for some reason
// Maybe related to https://github.com/rust-lang/rust/issues/47133
#![allow(dead_code)]  // Cheeky-breeky trick the linter

// Euler
use std::f64::consts::E;

// Constants used in definition of the PDFs
const INV_SQRT_2_PI: f64 = 0.3989422804014;   // (2 pi)^(-1/2)
const P1: [f64; 5] = [0.4259894875, -0.1249762550, 0.03984243700, -0.006298287635, 0.001511162253];
const Q1: [f64; 5] = [1.0, -0.3388260629, 0.09594393323, -0.01608042283, 0.003778942063];
const P2: [f64; 5] = [0.1788541609, 0.1173957403, 0.01488850518, -0.001394989411, 0.0001283617211];
const Q2: [f64; 5] = [1.0, 0.7428795082, 0.3153932961, 0.06694219548, 0.008790609714];
const P3: [f64; 5] = [0.1788544503, 0.09359161662, 0.006325387654, 0.00006611667319, -0.000002031049101];
const Q3: [f64; 5] = [1.0, 0.6097809921, 0.2560616665, 0.04746722384, 0.006957301675];
const P4: [f64; 5] = [0.9874054407, 118.6723273, 849.2794360, -743.7792444, 427.0262186];
const Q4: [f64; 5] = [1.0, 106.8615961, 337.6496214, 2016.712389, 1597.063511];
const P5: [f64; 5] = [1.003675074, 167.5702434, 4789.711289, 21217.86767, -22324.94910];
const Q5: [f64; 5] = [1.0, 156.9424537, 3745.310488, 9834.698876, 66924.28357];
const P6: [f64; 5] = [1.000827619, 664.9143136, 62972.92665, 475554.6998, -5743609.109];
const Q6: [f64; 5] = [1.0, 651.4101098, 56974.73333, 165917.4725, -2815759.939];
const A1: [f64; 3] = [0.04166666667, -0.01996527778, 0.02709538966];
const A2: [f64; 2] = [-1.845568670, -4.284640743];

/// Landau PDF according to GNU Scientific Library (GSL)
pub fn landau_pdf(x: f64, x0: f64, xi: f64) -> f64 {
    
    // Return immediately if our xi is le 0
    if xi <= 0.0 {
        return 0.0;
    }

    // Quantifier deciding the branch we go in
    let v = (x -x0) / xi;

    // Variables to be set
    let (u, ue, us, denlan): (f64, f64, f64, f64);

    if v < -5.5 {
        u = E.powf(v + 1.0);
        if u < 1e-10 {
            return 0.0;
        }
        ue = E.powf(-1.0 / u);
        us = u.sqrt();
        denlan = INV_SQRT_2_PI * (ue / us) * (1.0 + (A1[0] + (A1[1] + A1[2] * u)* u) *u);
    } else if  v < -1.0 {
        u = E.powf(-v -1.0);
        denlan = E.powf(-u) * u.sqrt() * (P1[0] + (P1[1] + (P1[2] + (P1[3] + P1[4] * v) *v) *v) *v) / (Q1[0] + (Q1[1] + (Q1[2] + (Q1[3] + Q1[4] * v) *v) *v) *v);
    } else if v < 1.0 {
        denlan = (P2[0] + (P2[1] + (P2[2] + (P2[3] + P2[4] * v) * v) * v) * v) / (Q2[0] + (Q2[1] + (Q2[2] + (Q2[3] + Q2[4] * v) * v) * v) * v);
    } else if v < 5.0 {
        denlan = (P3[0] + (P3[1] + (P3[2] + (P3[3] + P3[4] * v) * v) * v) * v) / (Q3[0] + (Q3[1] + (Q3[2] + (Q3[3] + Q3[4] * v) * v) * v) * v);
    } else if v < 12.0 {
        u = 1.0 / v;
		denlan = u * u * (P4[0] + (P4[1] + (P4[2] + (P4[3] + P4[4] * u) * u) * u) * u) / (Q4[0] + (Q4[1] + (Q4[2] + (Q4[3] + Q4[4] * u) * u) * u) * u);
    } else if v < 50.0 {
        u = 1.0 / v;
		denlan = u * u * (P5[0] + (P5[1] + (P5[2] + (P5[3] + P5[4] * u) * u) * u) * u) / (Q5[0] + (Q5[1] + (Q5[2] + (Q5[3] + Q5[4] * u) * u) * u) * u);
    } else if v < 300.0 {
        u = 1.0 / v;
		denlan = u * u * (P6[0] + (P6[1] + (P6[2] + (P6[3] + P6[4] * u) * u) * u) * u) / (Q6[0] + (Q6[1] + (Q6[2] + (Q6[3] + Q6[4] * u) * u) * u) * u);
    } else {
        u = 1.0 / (v - v * v.ln() / (v + 1.0));
		denlan = u * u * (1.0 + (A2[0] + A2[1] * u) * u);
    }
    // Return expression
    denlan / xi
}

/// Gaussian PDF
pub fn gauss_pdf(x: f64, mu: f64, sigma: f64) -> f64 {
    INV_SQRT_2_PI / sigma * E.powf(-0.5 * (x - mu).powf(2.0) / sigma.powf(2.0))
}

/// Convolution of Landau and Gaussian PDF
pub fn langau_pdf(x: f64, mu: f64, eta: f64, sigma: f64) -> f64 {

    // Control constants
    let mut n_steps = 100u32;  // Number of convolution steps
    let sigma_env = 8u8;  // Extend convolution to +- sigma_env Gaussian sigmas

    // Convolution steps have to be increased if sigma > eta * 5 to get stable solution that does not oscillate
    if sigma > 3.0 * eta {
        n_steps *= (sigma / eta / 3.0) as u32;
        if n_steps > 100000 {
            n_steps = 100000;}
    }

    // Range of convolution integral
	let x_low = x - sigma_env as f64 * sigma;
	let x_upp = x + sigma_env as f64 * sigma;
    let step = (x_upp - x_low) / n_steps as f64;

    // MP shift correction FIXME: This does nothing
    let mpshift: f64 = 0.0;  // -0.22278298;     // Landau maximum location shift in original code is wrong, since the shift does not depend on mu only
	let mpc = mu - mpshift;  // * eta;

    // Variables for putting in results
    let (mut xx, mut bound, mut sum): (f64, f64, f64);
    
    sum = 0.0;

    // Discrete linear convolution of Landau and Gaussian
    for i in 1..=(n_steps / 2) {
        bound = (i as f64 - 0.5) * step;
        xx = x_low + bound;
        sum += landau_pdf(xx, mpc, eta) / eta * gauss_pdf(x, xx, sigma);
        xx = x_upp - bound;
        sum += landau_pdf(xx, mpc, eta) / eta * gauss_pdf(x, xx, sigma);
    }
    
    // FIXME: Unsused in original implementation
    // const double norm = 0.398902310115109;  // normalization of the integral from -inf..intf
    
    // Return result
    step * sum
}