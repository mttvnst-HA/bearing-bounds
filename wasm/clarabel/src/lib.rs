//! Minimal WebAssembly binding of Clarabel for Bearing Bounds.
//!
//! Solves   minimise q'x   subject to   A x + s = b,  s in K
//! with K = zero^nz x nonneg^nl x SOC(d1) x SOC(d2) x ...  (no quadratic term).
//!
//! Unlike the cvxjs binding, the primal x, the dual z and the slack s are returned for every
//! status, including AlmostSolved. The caller decides what to accept: Bearing Bounds rebuilds
//! the field and checks it against the exact constraints before it reports a bound.

use clarabel::algebra::CscMatrix;
use clarabel::solver::*;
use wasm_bindgen::prelude::*;

/// Keep in step with the pinned version in Cargo.toml.
const CLARABEL_VERSION: &str = "0.11.1";

/// Status codes returned in out[0]. Keep in step with CLARABEL_STATUS in src/engine.js.
fn status_code(s: SolverStatus) -> f64 {
    match s {
        SolverStatus::Unsolved => 0.0,
        SolverStatus::Solved => 1.0,
        SolverStatus::PrimalInfeasible => 2.0,
        SolverStatus::DualInfeasible => 3.0,
        SolverStatus::AlmostSolved => 4.0,
        SolverStatus::AlmostPrimalInfeasible => 5.0,
        SolverStatus::AlmostDualInfeasible => 6.0,
        SolverStatus::MaxIterations => 7.0,
        SolverStatus::MaxTime => 8.0,
        SolverStatus::NumericalError => 9.0,
        SolverStatus::InsufficientProgress => 10.0,
        SolverStatus::CallbackTerminated => 11.0,
    }
}

#[wasm_bindgen]
pub fn version() -> String {
    format!("clarabel-bb {} (Clarabel.rs {})", env!("CARGO_PKG_VERSION"), CLARABEL_VERSION)
}

/// settings = [max_iter, time_limit_s, tol_gap_abs, tol_gap_rel, tol_feas, tol_ktratio, verbose,
///             equilibrate (0/1), presolve (0/1), static_regularization_constant,
///             dynamic_regularization_eps, dynamic_regularization_delta, iterative_refinement_reltol,
///             iterative_refinement_abstol, iterative_refinement_max_iter, max_step_fraction,
///             reduced_tol_gap_abs, reduced_tol_gap_rel, reduced_tol_feas, reduced_tol_ktratio]
/// (missing entries take Clarabel's defaults)
/// Returns [status, iterations, solve_time_s, obj_primal, obj_dual, x (n), z (m), s (m)].
/// On a setup error returns [-1] (bad data).
#[wasm_bindgen]
pub fn solve(
    n: usize,
    m: usize,
    q: &[f64],
    a_colptr: &[u32],
    a_rowidx: &[u32],
    a_val: &[f64],
    b: &[f64],
    n_zero: usize,
    n_nonneg: usize,
    soc_dims: &[u32],
    settings: &[f64],
) -> Vec<f64> {
    let p = CscMatrix::<f64>::zeros((n, n));
    let a = CscMatrix::new(
        m,
        n,
        a_colptr.iter().map(|&v| v as usize).collect(),
        a_rowidx.iter().map(|&v| v as usize).collect(),
        a_val.to_vec(),
    );
    let mut cones: Vec<SupportedConeT<f64>> = Vec::with_capacity(2 + soc_dims.len());
    if n_zero > 0 {
        cones.push(ZeroConeT(n_zero));
    }
    if n_nonneg > 0 {
        cones.push(NonnegativeConeT(n_nonneg));
    }
    for &d in soc_dims {
        cones.push(SecondOrderConeT(d as usize));
    }
    let g = |i: usize, d: f64| settings.get(i).copied().unwrap_or(d);
    let st = DefaultSettingsBuilder::<f64>::default()
        .max_iter(g(0, 200.0) as u32)
        .time_limit(g(1, f64::INFINITY))
        .tol_gap_abs(g(2, 1e-8))
        .tol_gap_rel(g(3, 1e-8))
        .tol_feas(g(4, 1e-8))
        .tol_ktratio(g(5, 1e-6))
        .verbose(g(6, 0.0) != 0.0)
        .equilibrate_enable(g(7, 1.0) != 0.0)
        .presolve_enable(g(8, 1.0) != 0.0)
        .static_regularization_constant(g(9, 1e-8))
        .dynamic_regularization_eps(g(10, 1e-13))
        .dynamic_regularization_delta(g(11, 2e-7))
        .iterative_refinement_reltol(g(12, 1e-13))
        .iterative_refinement_abstol(g(13, 1e-12))
        .iterative_refinement_max_iter(g(14, 10.0) as u32)
        .max_step_fraction(g(15, 0.99))
        .reduced_tol_gap_abs(g(16, 5e-5))
        .reduced_tol_gap_rel(g(17, 5e-5))
        .reduced_tol_feas(g(18, 1e-4))
        .reduced_tol_ktratio(g(19, 1e-4))
        .build()
        .unwrap();
    let mut solver = match DefaultSolver::new(&p, q, &a, b, &cones, st) {
        Ok(s) => s,
        Err(_) => return vec![-1.0],
    };
    solver.solve();
    let sol = &solver.solution;
    let mut out = Vec::with_capacity(5 + n + 2 * m);
    out.push(status_code(sol.status));
    out.push(sol.iterations as f64);
    out.push(sol.solve_time);
    out.push(sol.obj_val);
    out.push(sol.obj_val_dual);
    out.extend_from_slice(&sol.x);
    out.extend_from_slice(&sol.z);
    out.extend_from_slice(&sol.s);
    out
}
