/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@ssec.wisc.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 2.  See the file COPYING for more details.
**
*******************************************************************************/

%{
#include <gutil.h>

#include "input_local.h"

void comment(locus_data *locus);
%}


%option never-interactive
%option nounistd
%option noyywrap


%%
[ \t]+					;

"/*"					{ comment(locus); }
"//"[^\n]*              		{ /* consume //-comment */ }

-?fd					{ locus->statement = "fd";              return INPUT_OPTION_FD; }

-?options				{ locus->statement = "options";         return XRTM_INPUT_OPTIONS; }
-?solvers				{ locus->statement = "solvers";         return XRTM_INPUT_SOLVERS; }
-?max_coef				{ locus->statement = "max_coef";        return XRTM_INPUT_MAX_COEF; }
-?n_quad				{ locus->statement = "n_quad";          return XRTM_INPUT_N_QUAD; }
-?n_stokes				{ locus->statement = "n_stokes";        return XRTM_INPUT_N_STOKES; }
-?n_derivs				{ locus->statement = "n_derivs";        return XRTM_INPUT_N_DERIVS; }
-?n_layers				{ locus->statement = "n_layers";        return XRTM_INPUT_N_LAYERS; }
-?n_kernels				{ locus->statement = "n_kernels";       return XRTM_INPUT_N_KERNELS; }
-?n_kernel_quad				{ locus->statement = "n_kernel_quad";   return XRTM_INPUT_N_KERNEL_QUAD; }
-?kernels				{ locus->statement = "kernels";         return XRTM_INPUT_KERNELS; }
-?n_out_levels				{ locus->statement = "n_out_levels";    return XRTM_INPUT_N_OUT_LEVELS; }
-?n_out_thetas				{ locus->statement = "n_out_thetas";    return XRTM_INPUT_N_OUT_THETAS; }

-?doub_d_tau				{ locus->statement = "doub_d_tau";      return XRTM_INPUT_DOUB_D_TAU; }
-?pade_params				{ locus->statement = "pade_params";     return XRTM_INPUT_PADE_PARAMS; }
-?sos_params				{ locus->statement = "sos_params";      return XRTM_INPUT_SOS_PARAMS; }
-?fourier_tol				{ locus->statement = "fourier_tol";     return XRTM_INPUT_FOURIER_TOL; }
-?lambda				{ locus->statement = "lambda";          return XRTM_INPUT_LAMBDA; }
-?F_0					{ locus->statement = "F_0";             return XRTM_INPUT_F_0; }
-?theta_0				{ locus->statement = "theta_0";         return XRTM_INPUT_THETA_0; }
-?phi_0					{ locus->statement = "phi_0";           return XRTM_INPUT_PHI_0; }
-?out_levels				{ locus->statement = "out_levels";      return XRTM_INPUT_OUT_LEVELS; }
-?out_taus				{ locus->statement = "out_taus";        return XRTM_INPUT_OUT_TAUS; }
-?out_thetas				{ locus->statement = "out_thetas";      return XRTM_INPUT_OUT_THETAS; }
-?top_b					{ locus->statement = "top_b";           return XRTM_INPUT_TOP_B; }
-?top_b_l				{ locus->statement = "top_b_l";         return XRTM_INPUT_TOP_B_L; }
-?planet_r				{ locus->statement = "planet_r";        return XRTM_INPUT_PLANET_R; }
-?levels_z				{ locus->statement = "levels_z";        return XRTM_INPUT_LEVELS_Z; }
-?levels_b				{ locus->statement = "levels_b";        return XRTM_INPUT_LEVELS_B; }
-?levels_b_l				{ locus->statement = "levels_b_l";      return XRTM_INPUT_LEVELS_B_L; }
-?g					{ locus->statement = "g";               return XRTM_INPUT_G; }
-?g_l					{ locus->statement = "g_l";             return XRTM_INPUT_G_L; }
-?coef					{ locus->statement = "coef";            return XRTM_INPUT_CHI; }
-?coef_l				{ locus->statement = "coef_l";          return XRTM_INPUT_COEF_L; }
-?coef_files				{ locus->statement = "coef_files";      return XRTM_INPUT_COEF_FILES; }
-?coef_files_l				{ locus->statement = "coef_files_l";    return XRTM_INPUT_COEF_FILES_L; }
-?omega					{ locus->statement = "omega";           return XRTM_INPUT_OMEGA; }
-?omega_l				{ locus->statement = "omega_l";         return XRTM_INPUT_OMEGA_L; }
-?ltau					{ locus->statement = "ltau";            return XRTM_INPUT_LTAU; }
-?ltau_l				{ locus->statement = "ltau_l";          return XRTM_INPUT_LTAU_L; }
-?surface_b				{ locus->statement = "surface_b";       return XRTM_INPUT_SURFACE_B; }
-?surface_b_l				{ locus->statement = "surface_b-l";     return XRTM_INPUT_SURFACE_B_L; }
-?kernel_ampfac				{ locus->statement = "kernel_ampfac";   return XRTM_INPUT_KERNEL_AMPFAC; }
-?kernel_params				{ locus->statement = "kernel_params";   return XRTM_INPUT_KERNEL_PARAMS; }
-?kernel_ampfac_l			{ locus->statement = "kernel_ampfac_l"; return XRTM_INPUT_KERNEL_AMPFAC_L; }
-?kernel_params_l			{ locus->statement = "kernel_params_l"; return XRTM_INPUT_KERNEL_PARAMS_L; }

-?delta					{ locus->statement = "delta";           return XRTM_INPUT_DELTA; }
-?top_b_p				{ locus->statement = "top_b_p";         return XRTM_INPUT_TOP_B_P; }
-?levels_b_p				{ locus->statement = "levels_b_p";      return XRTM_INPUT_LEVELS_B_P; }
-?g_p					{ locus->statement = "g_p";             return XRTM_INPUT_G_P; }
-?coef_p				{ locus->statement = "coef_p";          return XRTM_INPUT_COEF_P; }
-?coef_files_p				{ locus->statement = "coef_files_p";    return XRTM_INPUT_COEF_FILES_P; }
-?omega_p				{ locus->statement = "omega_p";         return XRTM_INPUT_OMEGA_P; }
-?ltau_p				{ locus->statement = "ltau_p";          return XRTM_INPUT_LTAU_P; }
-?surface_b_p				{ locus->statement = "surface_b_p";     return XRTM_INPUT_SURFACE_B_P; }
-?kernel_ampfac_p			{ locus->statement = "kernel_ampfac_p"; return XRTM_INPUT_KERNEL_AMPFAC_P; }
-?kernel_params_p			{ locus->statement = "kernel_params_p"; return XRTM_INPUT_KERNEL_PARAMS_P; }
 /*
[A-Za-z_][-A-Za-z0-9_]*		{
     lex_type->s = strdup(yytext);
     return LEX_TYPE_IDENTIFIER;
}
 */
\'[^\'\n]*\' |
\"[^\"\n]*\"				{
     char c;
     c = yytext[yyleng - 1];
     yytext[yyleng - 1] = '\0';
     lex_type->s = strdup(yytext + 1);
     yytext[yyleng - 1] = c;
     return LEX_TYPE_STRING;
}

-?[0-9]+"."[0-9]+([eE][+-]?[0-9]+)?	|
-?[0-9]+"."([eE][+-]?[0-9]+)?		|
-?"."[0-9]+([eE][+-]?[0-9]+)?		{
     lex_type->d = atof(yytext);
     return LEX_TYPE_DOUBLE;
}

-?[0-9]+				{
     lex_type->i = atoi(yytext);
     return LEX_TYPE_INT;
}

[A-Za-z_][-A-Za-z0-9_:/.]*		{
     lex_type->s = strdup(yytext);
     return LEX_TYPE_STRING;
}

\n					{ locus->line++; }

.					{ return yytext[0]; }

%%


char *get_yytext() {

     return yytext;
}



void comment(locus_data *locus) {
	char c, prev = 0;
#ifdef __cplusplus
	while ((c = yyinput()) != 0)      /* (EOF maps to 0) */
#else
	while ((c =   input()) != 0)      /* (EOF maps to 0) */
#endif
	{
		if (c == '/' && prev == '*')
			return;
                else
                if (c == '\n')
                     locus->line++;
		prev = c;
	}
	input_error(locus, "unterminated comment");
}



void yypreinclude(FILE *fp) {

     yypush_buffer_state(yy_create_buffer(fp, YY_BUF_SIZE));
}



void yypostinclude() {

     yypop_buffer_state();
}



void yyrewind(int r, lex_type_data *type) {

    int i;

    for (i = yyleng - 1; i >= 0; --i)
        unput(yytext[i]);

    if (r == LEX_TYPE_STRING || r == LEX_TYPE_IDENTIFIER)
         free(type->s);
}
