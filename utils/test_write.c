/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>

#include <xrtm_interface.h>

#include "test.h"
#include "test_result.h"


/*******************************************************************************
 *
 ******************************************************************************/
static int write_result_header_bin(FILE *fp, int index, int index2, int i_out_level, int i_deriv, double out_theta, double out_phi, int i_stokes) {

     fwrite(&index,       sizeof(int),    1, fp);
     fwrite(&index2,      sizeof(int),    1, fp);
     fwrite(&i_out_level, sizeof(int),    1, fp);
     fwrite(&i_deriv,     sizeof(int),    1, fp);
     fwrite(&out_theta,   sizeof(double), 1, fp);
     fwrite(&out_phi,     sizeof(double), 1, fp);
     fwrite(&i_stokes,    sizeof(int),    1, fp);

     return 0;
}



static int write_result_header_text(FILE *fp, int index, int index2, int i_out_level, int i_deriv, double out_theta, double out_phi, int i_stokes) {

     if (i_deriv < 0)
          fprintf(fp, "%06d, %06d, i_out_level = %d,                out_theta = %9f, out_phi = %10f, i_stokes = %d: ", index, index2, i_out_level,          out_theta, out_phi, i_stokes);
     else
          fprintf(fp, "%06d, %06d, i_out_level = %d, i_deriv = %3d, out_theta = %9f, out_phi = %10f, i_stokes = %d: ", index, index2, i_out_level, i_deriv, out_theta, out_phi, i_stokes);

     return 0;
}



static int write_result_field_header_bin(FILE *fp, enum xrtm_solver_mask solver) {

     fwrite(&solver, sizeof(int), 1, fp);

     return 0;
}



static int write_result_field_header_blank_bin(FILE *fp, enum xrtm_solver_mask solver) {

     int minus_one = -1;

     fwrite(&minus_one, sizeof(int), 1, fp);

     return 0;
}



static int write_result_field_header_text(FILE *fp, enum xrtm_solver_mask solver) {

     fprintf(fp, "   %s -> ", xrtm_solver_mask_to_name(solver));

     return 0;
}



static int write_result_field_header_blank_text(FILE *fp, int i_solver) {

     int i;

     int length;

     char temp[1024];

     length = 3 + strlen(xrtm_solver_index_to_name(i_solver)) + 4;

     for (i = 0; i < length; ++i)
          temp[i] = ' ';

     temp[length] = '\0';
     fprintf(fp, "%s", temp);

     return 0;
}



static int write_result_field_values_bin(FILE *fp, double x_p, double x_m) {

     fwrite(&x_p, sizeof(double), 1, fp);
     fwrite(&x_m, sizeof(double), 1, fp);

     return 0;
}



static int write_result_field_values_blank_bin(FILE *fp) {

     return 0;
}



static int write_result_field_values_text(FILE *fp, double x_p, double x_m) {

     fprintf(fp, " % 13e % 13e", x_p, x_m);

     return 0;
}



static int write_result_field_values_blank_text(FILE *fp) {

     int i;

     int length;

     char temp[1024];

     length = 2 * 14;

     for (i = 0; i < length; ++i)
          temp[i] = ' ';

     temp[length] = '\0';
     fprintf(fp, "%s", temp);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int write_cmds(test_data *td) {

     test_cmd_data *test_cmd;

     test_xrtm_data *test_xrtm;

     list_for_each(&td->test_xrtm_list, test_xrtm) {
          list_for_each(&test_xrtm->test_cmd_list, test_cmd)
               fprintf(td->fp_core_xrtm_cmd, "%06d: %s\n", test_cmd->index, test_cmd->cmd);
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int write_results(test_data *td) {

     int i_solver;
     int i_solver2;
     int n_solvers;

     test_result_data *test_result;

     test_xrtm_data *test_xrtm;

     n_solvers = xrtm_solver_n();

     list_for_each(&td->test_xrtm_list, test_xrtm) {
          list_for_each(&test_xrtm->test_result_list, test_result) {
               if (td->core_write_results_bin) {
                    write_result_header_bin (td->fp_core_results_bin,           test_result->index,  test_result->index2, test_result->i_out_level, test_result->i_deriv, test_result->out_theta, test_result->out_phi, test_result->i_stokes);
                    write_result_header_bin (td->fp_core_results_bin_w_deltas,  test_result->index,  test_result->index2, test_result->i_out_level, test_result->i_deriv, test_result->out_theta, test_result->out_phi, test_result->i_stokes);
               }
               if (td->core_write_results_text) {
                    write_result_header_text(td->fp_core_results_text,          test_result->index,  test_result->index2, test_result->i_out_level, test_result->i_deriv, test_result->out_theta, test_result->out_phi, test_result->i_stokes);
                    write_result_header_text(td->fp_core_results_text_w_deltas, test_result->index,  test_result->index2, test_result->i_out_level, test_result->i_deriv, test_result->out_theta, test_result->out_phi, test_result->i_stokes);
               }

               for (i_solver = 0; i_solver < n_solvers; ++i_solver) {
                    for (i_solver2 = 0; i_solver2 < test_result->n_solvers; ++i_solver2) {
                         if (xrtm_solver_index_to_mask(i_solver) == test_result->solvers[i_solver2].mask)
                              break;
                    }

                    if (i_solver2 == test_result->n_solvers) {
                         if (td->core_fill_blanks) {
                              if (td->core_write_results_bin) {
                                   write_result_field_header_blank_bin(td->fp_core_results_bin, i_solver);
                                   write_result_field_values_blank_bin(td->fp_core_results_bin);

                                   write_result_field_header_blank_bin(td->fp_core_results_bin_w_deltas, i_solver);
                                   write_result_field_values_blank_bin(td->fp_core_results_bin_w_deltas);
                              }

                              if (td->core_write_results_text) {
                                   write_result_field_header_blank_text(td->fp_core_results_text, i_solver);
                                   write_result_field_values_blank_text(td->fp_core_results_text);

                                   write_result_field_header_blank_text(td->fp_core_results_text_w_deltas, i_solver);
                                   write_result_field_values_blank_text(td->fp_core_results_text_w_deltas);
                              }
                         }

                         continue;
                    }

                    if (td->core_write_results_bin) {
                         write_result_field_header_bin(td->fp_core_results_bin,   test_result->solvers[i_solver2].mask);
                         write_result_field_values_bin(td->fp_core_results_bin,   test_result->solvers[i_solver2].x_p, test_result->solvers[i_solver2].x_m);
                    }
                    if (td->core_write_results_text) {
                         write_result_field_header_text(td->fp_core_results_text, test_result->solvers[i_solver2].mask);
                         write_result_field_values_text(td->fp_core_results_text, test_result->solvers[i_solver2].x_p, test_result->solvers[i_solver2].x_m);
                    }

                    if (td->core_check_diffs) {
                         if (td->core_write_results_bin) {
                              write_result_field_header_bin(td->fp_core_results_bin_w_deltas, test_result->solvers[i_solver2].mask);
                              write_result_field_values_bin(td->fp_core_results_bin_w_deltas, test_result->solvers[i_solver2].x_p_delta, test_result->solvers[i_solver2].x_m_delta);
                         }
                         if (td->core_write_results_text) {
                              write_result_field_header_text(td->fp_core_results_text_w_deltas, test_result->solvers[i_solver2].mask);
                              write_result_field_values_text(td->fp_core_results_text_w_deltas, test_result->solvers[i_solver2].x_p_delta, test_result->solvers[i_solver2].x_m_delta);
                         }
                    }
               }

               if (td->core_write_results_text) {
                    fprintf(td->fp_core_results_text,          "\n");
                    fprintf(td->fp_core_results_text_w_deltas, "\n");
               }
          }
     }

     return 0;
}
