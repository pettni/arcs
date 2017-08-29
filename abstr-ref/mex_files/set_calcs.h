
#ifndef SET_CALCS
#define SET_CALCS

#define WIN_SET 1
#define WIN_CANDIDATE_SET 2
#define WIN_CANDIDATE_CONT 3

DdNode** win_primal(BDDSys* sys, DdNode* A, DdNode* B, DdNode** C, int C_num,
                   char quant, char quant2, DdNode* head_start, int mode);

DdNode** win_intermediate(BDDSys* sys, DdNode* A, DdNode* B, DdNode* Z,
                         DdNode** C, int C_num, char quant, int mode);
DdNode** win_until(BDDSys* sys, DdNode* Z, DdNode* B, char quant, int mode);
DdNode** win_until_and_always(BDDSys* sys, DdNode* A, DdNode* B, DdNode* Z,
                             char quant, int mode);
DdNode** PGpre(BDDSys* sys, DdNode* Z, DdNode* B, char quant, int mode);
DdNode** inv(BDDSys* sys, DdNode* Z, DdNode* B, DdNode* U, DdNode* G,
            char quant, int mode);

DdNode* pre(BDDSys* sys, DdNode* X_prime, DdNode* A, char quant1, char quant2);

#endif
