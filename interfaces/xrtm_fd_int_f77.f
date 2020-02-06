c***********************************************************************
c This file was generated by bindx version 0.01.  Edit at your own risk.
c***********************************************************************


      subroutine xrtm_fd_create_f77(d, xrtm, n_derivs, error)
      implicit none
      byte d(*)
      byte xrtm(*)
      integer n_derivs
      integer error
      integer xrtm_fd_create_bindx_f77
      error = xrtm_fd_create_bindx_f77(d, xrtm, n_derivs)
      end subroutine xrtm_fd_create_f77


      subroutine xrtm_fd_destroy_f77(d, error)
      implicit none
      byte d(*)
      integer error
      integer xrtm_fd_destroy_bindx_f77
      error = xrtm_fd_destroy_bindx_f77(d)
      end subroutine xrtm_fd_destroy_f77


      integer function xrtm_fd_get_n_derivs_f77(d)
      implicit none
      byte d(*)
      integer xrtm_fd_get_n_derivs_bindx_f77
      xrtm_fd_get_n_derivs_f77 = xrtm_fd_get_n_derivs_bindx_f77(d)
      end function xrtm_fd_get_n_derivs_f77


      subroutine xrtm_fd_set_delta_f77(d, delta, error)
      implicit none
      byte d(*)
      real*8 delta(*)
      integer error
      integer xrtm_fd_set_delta_bindx_f77
      error = xrtm_fd_set_delta_bindx_f77(d, delta)
      end subroutine xrtm_fd_set_delta_f77


      subroutine xrtm_fd_get_delta_f77(d, delta, error)
      implicit none
      byte d(*)
      real*8 delta(*)
      integer error
      integer xrtm_fd_get_delta_bindx_f77
      error = xrtm_fd_get_delta_bindx_f77(d, delta)
      end subroutine xrtm_fd_get_delta_f77


      subroutine xrtm_fd_set_F_iso_top_p_1_f77(d, i_deriv, F_iso_top_l,
     &error)
      implicit none
      byte d(*)
      integer i_deriv
      real*8 F_iso_top_l
      integer error
      integer xrtm_fd_set_F_iso_top_p_1_bindx_f77
      error = xrtm_fd_set_F_iso_top_p_1_bindx_f77(d, i_deriv,
     &F_iso_top_l)
      end subroutine xrtm_fd_set_F_iso_top_p_1_f77


      subroutine xrtm_fd_set_F_iso_top_p_n_f77(d, F_iso_top_l, error)
      implicit none
      byte d(*)
      real*8 F_iso_top_l(*)
      integer error
      integer xrtm_fd_set_F_iso_top_p_n_bindx_f77
      error = xrtm_fd_set_F_iso_top_p_n_bindx_f77(d, F_iso_top_l)
      end subroutine xrtm_fd_set_F_iso_top_p_n_f77


      real*8 function xrtm_fd_get_F_iso_top_p_f77(d, i_deriv)
      implicit none
      byte d(*)
      integer i_deriv
      integer xrtm_fd_get_F_iso_top_p_bindx_f77
      xrtm_fd_get_F_iso_top_p_f77 = xrtm_fd_get_F_iso_top_p_bindx_f77(d,
     &i_deriv)
      end function xrtm_fd_get_F_iso_top_p_f77


      subroutine xrtm_fd_set_F_iso_bot_p_1_f77(d, i_deriv, F_iso_bot_l,
     &error)
      implicit none
      byte d(*)
      integer i_deriv
      real*8 F_iso_bot_l
      integer error
      integer xrtm_fd_set_F_iso_bot_p_1_bindx_f77
      error = xrtm_fd_set_F_iso_bot_p_1_bindx_f77(d, i_deriv,
     &F_iso_bot_l)
      end subroutine xrtm_fd_set_F_iso_bot_p_1_f77


      subroutine xrtm_fd_set_F_iso_bot_p_n_f77(d, F_iso_bot_l, error)
      implicit none
      byte d(*)
      real*8 F_iso_bot_l(*)
      integer error
      integer xrtm_fd_set_F_iso_bot_p_n_bindx_f77
      error = xrtm_fd_set_F_iso_bot_p_n_bindx_f77(d, F_iso_bot_l)
      end subroutine xrtm_fd_set_F_iso_bot_p_n_f77


      real*8 function xrtm_fd_get_F_iso_bot_p_f77(d, i_deriv)
      implicit none
      byte d(*)
      integer i_deriv
      integer xrtm_fd_get_F_iso_bot_p_bindx_f77
      xrtm_fd_get_F_iso_bot_p_f77 = xrtm_fd_get_F_iso_bot_p_bindx_f77(d,
     &i_deriv)
      end function xrtm_fd_get_F_iso_bot_p_f77


      subroutine xrtm_fd_set_levels_b_p_1_f77(d, i_deriv, levels_b_l,
     &error)
      implicit none
      byte d(*)
      integer i_deriv
      real*8 levels_b_l(*)
      integer error
      integer xrtm_fd_set_levels_b_p_1_bindx_f77
      error = xrtm_fd_set_levels_b_p_1_bindx_f77(d, i_deriv, levels_b_l)
      end subroutine xrtm_fd_set_levels_b_p_1_f77


      subroutine xrtm_fd_set_levels_b_p_n_f77(d, levels_b_l, error)
      implicit none
      byte d(*)
      real*8 levels_b_l(*)
      integer error
      integer xrtm_fd_set_levels_b_p_n_bindx_f77
      error = xrtm_fd_set_levels_b_p_n_bindx_f77(d, levels_b_l)
      end subroutine xrtm_fd_set_levels_b_p_n_f77


      real*8 function xrtm_fd_get_levels_b_p_f77(d, i_deriv, levels_b_l)
      implicit none
      byte d(*)
      integer i_deriv
      real*8 levels_b_l(*)
      integer xrtm_fd_get_levels_b_p_bindx_f77
      xrtm_fd_get_levels_b_p_f77 = xrtm_fd_get_levels_b_p_bindx_f77(d,
     &i_deriv, levels_b_l)
      end function xrtm_fd_get_levels_b_p_f77


      subroutine xrtm_fd_set_surface_b_p_1_f77(d, i_deriv, surface_b_l,
     &error)
      implicit none
      byte d(*)
      integer i_deriv
      real*8 surface_b_l
      integer error
      integer xrtm_fd_set_surface_b_p_1_bindx_f77
      error = xrtm_fd_set_surface_b_p_1_bindx_f77(d, i_deriv,
     &surface_b_l)
      end subroutine xrtm_fd_set_surface_b_p_1_f77


      subroutine xrtm_fd_set_surface_b_p_n_f77(d, surface_b_l, error)
      implicit none
      byte d(*)
      real*8 surface_b_l(*)
      integer error
      integer xrtm_fd_set_surface_b_p_n_bindx_f77
      error = xrtm_fd_set_surface_b_p_n_bindx_f77(d, surface_b_l)
      end subroutine xrtm_fd_set_surface_b_p_n_f77


      real*8 function xrtm_fd_get_surface_b_p_f77(d, i_deriv)
      implicit none
      byte d(*)
      integer i_deriv
      integer xrtm_fd_get_surface_b_p_bindx_f77
      xrtm_fd_get_surface_b_p_f77 = xrtm_fd_get_surface_b_p_bindx_f77(d,
     &i_deriv)
      end function xrtm_fd_get_surface_b_p_f77


      subroutine xrtm_fd_set_g_p_11_f77(d, i_layer, i_deriv, g_l, error)
      implicit none
      byte d(*)
      integer i_layer
      integer i_deriv
      real*8 g_l
      integer error
      integer xrtm_fd_set_g_p_11_bindx_f77
      error = xrtm_fd_set_g_p_11_bindx_f77(d, i_layer, i_deriv, g_l)
      end subroutine xrtm_fd_set_g_p_11_f77


      subroutine xrtm_fd_set_g_p_n1_f77(d, i_deriv, g_l, error)
      implicit none
      byte d(*)
      integer i_deriv
      real*8 g_l(*)
      integer error
      integer xrtm_fd_set_g_p_n1_bindx_f77
      error = xrtm_fd_set_g_p_n1_bindx_f77(d, i_deriv, g_l)
      end subroutine xrtm_fd_set_g_p_n1_f77


      subroutine xrtm_fd_set_g_p_1n_f77(d, i_layer, g_l, error)
      implicit none
      byte d(*)
      integer i_layer
      real*8 g_l(*)
      integer error
      integer xrtm_fd_set_g_p_1n_bindx_f77
      error = xrtm_fd_set_g_p_1n_bindx_f77(d, i_layer, g_l)
      end subroutine xrtm_fd_set_g_p_1n_f77


      subroutine xrtm_fd_set_g_p_nn_f77(d, g_l, error)
      implicit none
      byte d(*)
      real*8 g_l(*)
      integer error
      integer xrtm_fd_set_g_p_nn_bindx_f77
      error = xrtm_fd_set_g_p_nn_bindx_f77(d, g_l)
      end subroutine xrtm_fd_set_g_p_nn_f77


      real*8 function xrtm_fd_get_g_p_f77(d, i_layer, i_deriv)
      implicit none
      byte d(*)
      integer i_layer
      integer i_deriv
      integer xrtm_fd_get_g_p_bindx_f77
      xrtm_fd_get_g_p_f77 = xrtm_fd_get_g_p_bindx_f77(d, i_layer,
     &i_deriv)
      end function xrtm_fd_get_g_p_f77


      subroutine xrtm_fd_set_coef_p_11_f77(d, i_layer, i_deriv, coef_l,
     &error)
      implicit none
      byte d(*)
      integer i_layer
      integer i_deriv
      real*8 coef_l(*)
      integer error
      integer xrtm_fd_set_coef_p_11_bindx_f77
      error = xrtm_fd_set_coef_p_11_bindx_f77(d, i_layer, i_deriv,
     &coef_l)
      end subroutine xrtm_fd_set_coef_p_11_f77


      subroutine xrtm_fd_set_coef_p_n1_f77(d, i_deriv, coef_l, error)
      implicit none
      byte d(*)
      integer i_deriv
      real*8 coef_l(*)
      integer error
      integer xrtm_fd_set_coef_p_n1_bindx_f77
      error = xrtm_fd_set_coef_p_n1_bindx_f77(d, i_deriv, coef_l)
      end subroutine xrtm_fd_set_coef_p_n1_f77


      subroutine xrtm_fd_set_coef_p_1n_f77(d, i_layer, coef_l, error)
      implicit none
      byte d(*)
      integer i_layer
      real*8 coef_l(*)
      integer error
      integer xrtm_fd_set_coef_p_1n_bindx_f77
      error = xrtm_fd_set_coef_p_1n_bindx_f77(d, i_layer, coef_l)
      end subroutine xrtm_fd_set_coef_p_1n_f77


      subroutine xrtm_fd_set_coef_p_nn_f77(d, coef_l, error)
      implicit none
      byte d(*)
      real*8 coef_l(*)
      integer error
      integer xrtm_fd_set_coef_p_nn_bindx_f77
      error = xrtm_fd_set_coef_p_nn_bindx_f77(d, coef_l)
      end subroutine xrtm_fd_set_coef_p_nn_f77


      real*8 function xrtm_fd_get_coef_p_f77(d, i_layer, i_deriv,
     &i_elem, i_coef)
      implicit none
      byte d(*)
      integer i_layer
      integer i_deriv
      integer i_elem
      integer i_coef
      integer xrtm_fd_get_coef_p_bindx_f77
      xrtm_fd_get_coef_p_f77 = xrtm_fd_get_coef_p_bindx_f77(d, i_layer,
     &i_deriv, i_elem, i_coef)
      end function xrtm_fd_get_coef_p_f77


      subroutine xrtm_fd_set_omega_p_11_f77(d, i_layer, i_deriv,
     &omega_l, error)
      implicit none
      byte d(*)
      integer i_layer
      integer i_deriv
      real*8 omega_l
      integer error
      integer xrtm_fd_set_omega_p_11_bindx_f77
      error = xrtm_fd_set_omega_p_11_bindx_f77(d, i_layer, i_deriv,
     &omega_l)
      end subroutine xrtm_fd_set_omega_p_11_f77


      subroutine xrtm_fd_set_omega_p_n1_f77(d, i_deriv, omega_l, error)
      implicit none
      byte d(*)
      integer i_deriv
      real*8 omega_l(*)
      integer error
      integer xrtm_fd_set_omega_p_n1_bindx_f77
      error = xrtm_fd_set_omega_p_n1_bindx_f77(d, i_deriv, omega_l)
      end subroutine xrtm_fd_set_omega_p_n1_f77


      subroutine xrtm_fd_set_omega_p_1n_f77(d, i_layer, omega_l, error)
      implicit none
      byte d(*)
      integer i_layer
      real*8 omega_l(*)
      integer error
      integer xrtm_fd_set_omega_p_1n_bindx_f77
      error = xrtm_fd_set_omega_p_1n_bindx_f77(d, i_layer, omega_l)
      end subroutine xrtm_fd_set_omega_p_1n_f77


      subroutine xrtm_fd_set_omega_p_nn_f77(d, omega_l, error)
      implicit none
      byte d(*)
      real*8 omega_l(*)
      integer error
      integer xrtm_fd_set_omega_p_nn_bindx_f77
      error = xrtm_fd_set_omega_p_nn_bindx_f77(d, omega_l)
      end subroutine xrtm_fd_set_omega_p_nn_f77


      real*8 function xrtm_fd_get_omega_p_f77(d, i_layer, i_deriv)
      implicit none
      byte d(*)
      integer i_layer
      integer i_deriv
      integer xrtm_fd_get_omega_p_bindx_f77
      xrtm_fd_get_omega_p_f77 = xrtm_fd_get_omega_p_bindx_f77(d,
     &i_layer, i_deriv)
      end function xrtm_fd_get_omega_p_f77


      subroutine xrtm_fd_set_ltau_p_11_f77(d, i_layer, i_deriv, ltau_l,
     &error)
      implicit none
      byte d(*)
      integer i_layer
      integer i_deriv
      real*8 ltau_l
      integer error
      integer xrtm_fd_set_ltau_p_11_bindx_f77
      error = xrtm_fd_set_ltau_p_11_bindx_f77(d, i_layer, i_deriv,
     &ltau_l)
      end subroutine xrtm_fd_set_ltau_p_11_f77


      subroutine xrtm_fd_set_ltau_p_n1_f77(d, i_deriv, ltau_l, error)
      implicit none
      byte d(*)
      integer i_deriv
      real*8 ltau_l(*)
      integer error
      integer xrtm_fd_set_ltau_p_n1_bindx_f77
      error = xrtm_fd_set_ltau_p_n1_bindx_f77(d, i_deriv, ltau_l)
      end subroutine xrtm_fd_set_ltau_p_n1_f77


      subroutine xrtm_fd_set_ltau_p_1n_f77(d, i_layer, ltau_l, error)
      implicit none
      byte d(*)
      integer i_layer
      real*8 ltau_l(*)
      integer error
      integer xrtm_fd_set_ltau_p_1n_bindx_f77
      error = xrtm_fd_set_ltau_p_1n_bindx_f77(d, i_layer, ltau_l)
      end subroutine xrtm_fd_set_ltau_p_1n_f77


      subroutine xrtm_fd_set_ltau_p_nn_f77(d, ltau_l, error)
      implicit none
      byte d(*)
      real*8 ltau_l(*)
      integer error
      integer xrtm_fd_set_ltau_p_nn_bindx_f77
      error = xrtm_fd_set_ltau_p_nn_bindx_f77(d, ltau_l)
      end subroutine xrtm_fd_set_ltau_p_nn_f77


      real*8 function xrtm_fd_get_ltau_p_f77(d, i_layer, i_deriv)
      implicit none
      byte d(*)
      integer i_layer
      integer i_deriv
      integer xrtm_fd_get_ltau_p_bindx_f77
      xrtm_fd_get_ltau_p_f77 = xrtm_fd_get_ltau_p_bindx_f77(d, i_layer,
     &i_deriv)
      end function xrtm_fd_get_ltau_p_f77


      subroutine xrtm_fd_set_kernel_ampfac_p_1_f77(d, i_kernel, i_deriv,
     &ampfac_l, error)
      implicit none
      byte d(*)
      integer i_kernel
      integer i_deriv
      real*8 ampfac_l
      integer error
      integer xrtm_fd_set_kernel_ampfac_p_1_bindx_f77
      error = xrtm_fd_set_kernel_ampfac_p_1_bindx_f77(d, i_kernel,
     &i_deriv, ampfac_l)
      end subroutine xrtm_fd_set_kernel_ampfac_p_1_f77


      subroutine xrtm_fd_set_kernel_ampfac_p_n_f77(d, i_kernel,
     &ampfac_l, error)
      implicit none
      byte d(*)
      integer i_kernel
      real*8 ampfac_l(*)
      integer error
      integer xrtm_fd_set_kernel_ampfac_p_n_bindx_f77
      error = xrtm_fd_set_kernel_ampfac_p_n_bindx_f77(d, i_kernel,
     &ampfac_l)
      end subroutine xrtm_fd_set_kernel_ampfac_p_n_f77


      real*8 function xrtm_fd_get_kernel_ampfac_p_f77(d, i_kernel,
     &i_deriv)
      implicit none
      byte d(*)
      integer i_kernel
      integer i_deriv
      integer xrtm_fd_get_kernel_ampfac_p_bindx_f77
      xrtm_fd_get_kernel_ampfac_p_f77 =
     &xrtm_fd_get_kernel_ampfac_p_bindx_f77(d, i_kernel, i_deriv)
      end function xrtm_fd_get_kernel_ampfac_p_f77


      subroutine xrtm_fd_set_kernel_params_p_11_f77(d, i_kernel,
     &i_deriv, i_param, param_l, error)
      implicit none
      byte d(*)
      integer i_kernel
      integer i_deriv
      integer i_param
      real*8 param_l
      integer error
      integer xrtm_fd_set_kernel_params_p_11_bindx_f77
      error = xrtm_fd_set_kernel_params_p_11_bindx_f77(d, i_kernel,
     &i_deriv, i_param, param_l)
      end subroutine xrtm_fd_set_kernel_params_p_11_f77


      subroutine xrtm_fd_set_kernel_params_p_1n_f77(d, i_kernel,
     &i_deriv, params_l, error)
      implicit none
      byte d(*)
      integer i_kernel
      integer i_deriv
      real*8 params_l(*)
      integer error
      integer xrtm_fd_set_kernel_params_p_1n_bindx_f77
      error = xrtm_fd_set_kernel_params_p_1n_bindx_f77(d, i_kernel,
     &i_deriv, params_l)
      end subroutine xrtm_fd_set_kernel_params_p_1n_f77


      subroutine xrtm_fd_set_kernel_params_p_n1_f77(d, i_kernel,
     &i_param, params_l, error)
      implicit none
      byte d(*)
      integer i_kernel
      integer i_param
      real*8 params_l(*)
      integer error
      integer xrtm_fd_set_kernel_params_p_n1_bindx_f77
      error = xrtm_fd_set_kernel_params_p_n1_bindx_f77(d, i_kernel,
     &i_param, params_l)
      end subroutine xrtm_fd_set_kernel_params_p_n1_f77


      subroutine xrtm_fd_set_kernel_params_p_nn_f77(d, i_kernel,
     &params_l, error)
      implicit none
      byte d(*)
      integer i_kernel
      real*8 params_l(*)
      integer error
      integer xrtm_fd_set_kernel_params_p_nn_bindx_f77
      error = xrtm_fd_set_kernel_params_p_nn_bindx_f77(d, i_kernel,
     &params_l)
      end subroutine xrtm_fd_set_kernel_params_p_nn_f77


      real*8 function xrtm_fd_get_kernel_params_p_f77(d, i_kernel,
     &i_deriv, i_param)
      implicit none
      byte d(*)
      integer i_kernel
      integer i_deriv
      integer i_param
      integer xrtm_fd_get_kernel_params_p_bindx_f77
      xrtm_fd_get_kernel_params_p_f77 =
     &xrtm_fd_get_kernel_params_p_bindx_f77(d, i_kernel, i_deriv, i_param)
      end function xrtm_fd_get_kernel_params_p_f77


      subroutine xrtm_fd_solution_f77(d, solver, solutions, method,
     &n_out_phis, out_phis, I_p, I_m, K_p, K_m, mean_p, mean_m,
     &mean_p_l, mean_m_l, flux_p, flux_m, flux_p_l, flux_m_l, flux_div, flux_div_l, error)
      implicit none
      byte d(*)
      integer solver
      integer solutions
      integer method
      integer n_out_phis
      real*8 out_phis(*)
      real*8 I_p(*)
      real*8 I_m(*)
      real*8 K_p(*)
      real*8 K_m(*)
      real*8 mean_p(*)
      real*8 mean_m(*)
      real*8 mean_p_l(*)
      real*8 mean_m_l(*)
      real*8 flux_p(*)
      real*8 flux_m(*)
      real*8 flux_p_l(*)
      real*8 flux_m_l(*)
      real*8 flux_div(*)
      real*8 flux_div_l(*)
      integer error
      integer xrtm_fd_solution_bindx_f77
      error = xrtm_fd_solution_bindx_f77(d, solver, solutions, method,
     &n_out_phis, out_phis, I_p, I_m, K_p, K_m, mean_p, mean_m,
     &mean_p_l, mean_m_l, flux_p, flux_m, flux_p_l, flux_m_l, flux_div, flux_div_l)
      end subroutine xrtm_fd_solution_f77


      subroutine xrtm_fd_radiance_f77(d, solver, method, n_out_phis,
     &out_phis, I_p, I_m, K_p, K_m, error)
      implicit none
      byte d(*)
      integer solver
      integer method
      integer n_out_phis
      real*8 out_phis(*)
      real*8 I_p(*)
      real*8 I_m(*)
      real*8 K_p(*)
      real*8 K_m(*)
      integer error
      integer xrtm_fd_radiance_bindx_f77
      error = xrtm_fd_radiance_bindx_f77(d, solver, method, n_out_phis,
     &out_phis, I_p, I_m, K_p, K_m)
      end subroutine xrtm_fd_radiance_f77

