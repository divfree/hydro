#pragma once

////////////////////////////////
// NAMES OF PARAMETERS
////////////////////////////////


// bool parameters
#define _PARDISO "PARDISO"
#define _implicit_boundary_conditions "implicit_boundary_conditions"
#define _FDMV "FDMV"
#define _coefficients_from_n "coefficients_from_n"
#define _replace_first_equation "replace_first_equation"
#define _SIMPLE "SIMPLE"
#define _MMIM "MMIM"
#define _sample_comparison "sample_comparison"
#define _pressure_corr_enable "pressure_corr_enable"
#define _iter_history_mesh "iter_history_mesh"
#define _outflow_zero_dpdn "outflow_zero_dpdn"
#define _stat_s_enable "stat_s_enable"

// double parameters
#define _A_decomp_eps "A_decomp_eps"
#define _eigen_corr_p "eigen_corr_p"
#define _eigen_corr_vel "eigen_corr_vel"
#define _alpha_u "alpha_u"
#define _alpha_v "alpha_v"
#define _alpha_p "alpha_p"
#define _eps_upwind "eps_upwind"
#define _imp_phi "imp_phi"
#define _imp_theta "imp_theta"
#define _exp_phi "exp_phi"
#define _exp_theta "exp_theta"
#define _SA_threshold "SA_threshold"

// string parameters
#define _scheme "scheme"
#define _linear_solver_vel "linear_solver_vel"
#define _linear_solver_p "linear_solver_p"
#define _linear_solver "linear_solver"
#define _bound_cond_vel "bound_cond_vel"
#define _bound_cond_p "bound_cond_p"
#define _bound_cond "bound_cond"
#define _linear_solver_nop "linear_solver_nop"
#define _SIMPLE_MOD "SIMPLE_MOD"
#define _problem "problem"
#define _plt_title "plt_title"
#define _output_dir "output_dir"
#define _exp_name "name"

// grid_double parameters
#define _grid_sample_u "grid_sample_u"
#define _grid_sample_v "grid_sample_v"
#define _grid_sample_p "grid_sample_p"

// OTX values names
#define _otx_volume "0volume"
#define _otx_x "0x"
#define _otx_y "0y"
#define _otx_u "1u"
#define _otx_v "1v"
#define _otx_w "1w"
#define _otx_p "2p"
#define _otx_beta "3beta"
#define _otx_div "4div"
#define _otx_divf "4divf"
#define _otx_markers_density "5markers_density"


// OT values names
#define _ot_t "0t"
#define _ot_s "1s"
#define _ot_last_s "1_ot_last_s"
#define _ot_gamma "2gamma"
#define _ot_div "2div"
#define _ot_div_mean "2div_mean"
#define _ot_Rn_u "3Rn_u"
#define _ot_Rn_v "3Rn_v"
#define _ot_Rn_w "3Rn_w"
#define _ot_Rn_p "3Rn_p"
#define _ot_Rn_max "3Rn_max"
#define _ot_Rs_u "4Rs_u"
#define _ot_Rs_v "4Rs_v"
#define _ot_Rs_w "4Rs_w"
#define _ot_Rs_p "4Rs_p"
#define _ot_Rs_max "4Rs_max"
#define _ot_Rsample_u "5Rsample_u"
#define _ot_Rsample_v "5Rsample_v"
#define _ot_Rsample_p "5Rsample_p"
#define _ot_Rsample_max "5Rsample_max"
#define _ot_Cd "6Cd"
#define _ot_Cl "6Cl"
#define _ot_recirc "7recirc"
#define _ot_interval "81interval"
#define _ot_value "82value"
