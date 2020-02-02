/*******************************************************************************
**
**    Copyright (C) 1998-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>

#include <xrtm_fd_interface.h>
#include <xrtm_interface.h>
#include <xrtm_support.h>


int main() {

     char temp[LN];
     char format[NM];

     int i;
     int n;

     int length;

     FILE *fp;

     fp = fopen("xrtm_enumeration.def", "w");


     length = xrtm_option_max_name_length();
     length = MAX(length, xrtm_solver_max_name_length());
     length = MAX(length, xrtm_output_max_name_length());
     length = MAX(length, xrtm_kernel_max_name_length());
     length = MAX(length, xrtm_fd_method_max_name_length());


     fprintf(fp, "enumeration option_mask,\n");
     n = xrtm_option_n();

     for (i = 0; i < n; ++i) {
          strtoupper(xrtm_option_index_to_name(i), temp);
          snprintf(format, NM, "     OPTION_%%-%ds = %%d", length);
          fprintf(fp, format, temp, xrtm_option_index_to_mask(i));
          if (i < n - 1)
               fprintf(fp, ",\n");
     }

     fprintf(fp, ";\n");
     fprintf(fp, "\n");


     fprintf(fp, "enumeration solver_mask,\n");

     n = xrtm_solver_n();

     for (i = 0; i < n; ++i) {
          strtoupper(xrtm_solver_index_to_name(i), temp);
          snprintf(format, NM, "     SOLVER_%%-%ds = %%d", length);
          fprintf(fp, format, temp, xrtm_solver_index_to_mask(i));
          if (i < n - 1)
               fprintf(fp, ",\n");
     }

     fprintf(fp, ";\n");
     fprintf(fp, "\n");


     fprintf(fp, "enumeration output_mask,\n");

     n = xrtm_output_n();

     for (i = 0; i < n; ++i) {
          strtoupper(xrtm_output_index_to_name(i), temp);
          snprintf(format, NM, "     OUTPUT_%%-%ds = %%d", length);
          fprintf(fp, format, temp, xrtm_output_index_to_mask(i));
          if (i < n - 1)
               fprintf(fp, ",\n");
     }

     fprintf(fp, ";\n");
     fprintf(fp, "\n");


     fprintf(fp, "enumeration kernel_type,\n");

     n = xrtm_kernel_n();

     for (i = 0; i < n; ++i) {
          strtoupper(xrtm_kernel_index_to_name(i), temp);
          snprintf(format, NM, "     KERNEL_%%-%ds = %%d", length);
          fprintf(fp, format, temp, xrtm_kernel_index_to_value(i));
          if (i < n - 1)
               fprintf(fp, ",\n");
     }

     fprintf(fp, ";\n");
     fprintf(fp, "\n");


     fp = fopen("xrtm_structure.def", "w");

     fprintf(fp, "structure xrtm_type %lu;\n", sizeof(xrtm_data));

     fclose(fp);


     fp = fopen("xrtm_fd_enumeration.def", "w");

     fprintf(fp, "enumeration method_type,\n");

     n = xrtm_fd_method_n();

     for (i = 0; i < n; ++i) {
          strtoupper(xrtm_fd_method_index_to_name(i), temp);
          snprintf(format, NM, "     METHOD_%%-%ds = %%d", length);
          fprintf(fp, format, temp, xrtm_fd_method_index_to_value(i));
          if (i < n - 1)
               fprintf(fp, ",\n");
     }

     fprintf(fp, ";\n");
     fprintf(fp, "\n");

     fclose(fp);


     fp = fopen("xrtm_fd_structure.def", "w");

     fprintf(fp, "structure xrtm_fd_type %lu;\n", sizeof(xrtm_fd_data));

     fclose(fp);

     exit(0);
}
