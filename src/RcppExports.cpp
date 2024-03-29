// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// dgCMatToJMat
void dgCMatToJMat(Rcpp::S4 q, std::string fname, std::string mtype, std::string ctype, std::string valuetype, bool transpose, std::string comment);
RcppExport SEXP _scellpam_dgCMatToJMat(SEXP qSEXP, SEXP fnameSEXP, SEXP mtypeSEXP, SEXP ctypeSEXP, SEXP valuetypeSEXP, SEXP transposeSEXP, SEXP commentSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type q(qSEXP);
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type mtype(mtypeSEXP);
    Rcpp::traits::input_parameter< std::string >::type ctype(ctypeSEXP);
    Rcpp::traits::input_parameter< std::string >::type valuetype(valuetypeSEXP);
    Rcpp::traits::input_parameter< bool >::type transpose(transposeSEXP);
    Rcpp::traits::input_parameter< std::string >::type comment(commentSEXP);
    dgCMatToJMat(q, fname, mtype, ctype, valuetype, transpose, comment);
    return R_NilValue;
END_RCPP
}
// GetSeuratGroups
Rcpp::IntegerVector GetSeuratGroups(Rcpp::S4 q);
RcppExport SEXP _scellpam_GetSeuratGroups(SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(GetSeuratGroups(q));
    return rcpp_result_gen;
END_RCPP
}
// SceToJMat
void SceToJMat(Rcpp::NumericMatrix& M, std::string fname, Rcpp::Nullable<Rcpp::StringVector> rownames, Rcpp::Nullable<Rcpp::StringVector> colnames, std::string mtype, std::string ctype, std::string valuetype, bool transpose, std::string comment);
RcppExport SEXP _scellpam_SceToJMat(SEXP MSEXP, SEXP fnameSEXP, SEXP rownamesSEXP, SEXP colnamesSEXP, SEXP mtypeSEXP, SEXP ctypeSEXP, SEXP valuetypeSEXP, SEXP transposeSEXP, SEXP commentSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type M(MSEXP);
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::StringVector> >::type rownames(rownamesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::StringVector> >::type colnames(colnamesSEXP);
    Rcpp::traits::input_parameter< std::string >::type mtype(mtypeSEXP);
    Rcpp::traits::input_parameter< std::string >::type ctype(ctypeSEXP);
    Rcpp::traits::input_parameter< std::string >::type valuetype(valuetypeSEXP);
    Rcpp::traits::input_parameter< bool >::type transpose(transposeSEXP);
    Rcpp::traits::input_parameter< std::string >::type comment(commentSEXP);
    SceToJMat(M, fname, rownames, colnames, mtype, ctype, valuetype, transpose, comment);
    return R_NilValue;
END_RCPP
}
// BuildAbundanceMatrix
Rcpp::NumericMatrix BuildAbundanceMatrix(Rcpp::NumericVector clasif, Rcpp::IntegerVector gr, int expgroups);
RcppExport SEXP _scellpam_BuildAbundanceMatrix(SEXP clasifSEXP, SEXP grSEXP, SEXP expgroupsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type clasif(clasifSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type gr(grSEXP);
    Rcpp::traits::input_parameter< int >::type expgroups(expgroupsSEXP);
    rcpp_result_gen = Rcpp::wrap(BuildAbundanceMatrix(clasif, gr, expgroups));
    return rcpp_result_gen;
END_RCPP
}
// CsvToJMat
void CsvToJMat(std::string ifname, std::string ofname, std::string mtype, char csep, std::string ctype, std::string valuetype, bool transpose, std::string comment);
RcppExport SEXP _scellpam_CsvToJMat(SEXP ifnameSEXP, SEXP ofnameSEXP, SEXP mtypeSEXP, SEXP csepSEXP, SEXP ctypeSEXP, SEXP valuetypeSEXP, SEXP transposeSEXP, SEXP commentSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type ifname(ifnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type ofname(ofnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type mtype(mtypeSEXP);
    Rcpp::traits::input_parameter< char >::type csep(csepSEXP);
    Rcpp::traits::input_parameter< std::string >::type ctype(ctypeSEXP);
    Rcpp::traits::input_parameter< std::string >::type valuetype(valuetypeSEXP);
    Rcpp::traits::input_parameter< bool >::type transpose(transposeSEXP);
    Rcpp::traits::input_parameter< std::string >::type comment(commentSEXP);
    CsvToJMat(ifname, ofname, mtype, csep, ctype, valuetype, transpose, comment);
    return R_NilValue;
END_RCPP
}
// JMatToCsv
void JMatToCsv(std::string ifile, std::string csvfile, char csep, bool withquotes);
RcppExport SEXP _scellpam_JMatToCsv(SEXP ifileSEXP, SEXP csvfileSEXP, SEXP csepSEXP, SEXP withquotesSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type ifile(ifileSEXP);
    Rcpp::traits::input_parameter< std::string >::type csvfile(csvfileSEXP);
    Rcpp::traits::input_parameter< char >::type csep(csepSEXP);
    Rcpp::traits::input_parameter< bool >::type withquotes(withquotesSEXP);
    JMatToCsv(ifile, csvfile, csep, withquotes);
    return R_NilValue;
END_RCPP
}
// ScellpamSetDebug
void ScellpamSetDebug(bool deb, bool debparpam, bool debjmat);
RcppExport SEXP _scellpam_ScellpamSetDebug(SEXP debSEXP, SEXP debparpamSEXP, SEXP debjmatSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< bool >::type deb(debSEXP);
    Rcpp::traits::input_parameter< bool >::type debparpam(debparpamSEXP);
    Rcpp::traits::input_parameter< bool >::type debjmat(debjmatSEXP);
    ScellpamSetDebug(deb, debparpam, debjmat);
    return R_NilValue;
END_RCPP
}
// CalcAndWriteDissimilarityMatrix
void CalcAndWriteDissimilarityMatrix(std::string ifname, std::string ofname, std::string distype, std::string restype, std::string comment, int nthreads);
RcppExport SEXP _scellpam_CalcAndWriteDissimilarityMatrix(SEXP ifnameSEXP, SEXP ofnameSEXP, SEXP distypeSEXP, SEXP restypeSEXP, SEXP commentSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type ifname(ifnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type ofname(ofnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type distype(distypeSEXP);
    Rcpp::traits::input_parameter< std::string >::type restype(restypeSEXP);
    Rcpp::traits::input_parameter< std::string >::type comment(commentSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    CalcAndWriteDissimilarityMatrix(ifname, ofname, distype, restype, comment, nthreads);
    return R_NilValue;
END_RCPP
}
// FilterJMatByName
void FilterJMatByName(std::string fname, Rcpp::StringVector Gn, std::string filname, std::string namesat);
RcppExport SEXP _scellpam_FilterJMatByName(SEXP fnameSEXP, SEXP GnSEXP, SEXP filnameSEXP, SEXP namesatSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type Gn(GnSEXP);
    Rcpp::traits::input_parameter< std::string >::type filname(filnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type namesat(namesatSEXP);
    FilterJMatByName(fname, Gn, filname, namesat);
    return R_NilValue;
END_RCPP
}
// FilterBySilhouetteQuantile
Rcpp::List FilterBySilhouetteQuantile(Rcpp::NumericVector s, Rcpp::List L, std::string fallcounts, std::string ffilcounts, std::string falldissim, std::string ffildissim, float q, bool addcom);
RcppExport SEXP _scellpam_FilterBySilhouetteQuantile(SEXP sSEXP, SEXP LSEXP, SEXP fallcountsSEXP, SEXP ffilcountsSEXP, SEXP falldissimSEXP, SEXP ffildissimSEXP, SEXP qSEXP, SEXP addcomSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type L(LSEXP);
    Rcpp::traits::input_parameter< std::string >::type fallcounts(fallcountsSEXP);
    Rcpp::traits::input_parameter< std::string >::type ffilcounts(ffilcountsSEXP);
    Rcpp::traits::input_parameter< std::string >::type falldissim(falldissimSEXP);
    Rcpp::traits::input_parameter< std::string >::type ffildissim(ffildissimSEXP);
    Rcpp::traits::input_parameter< float >::type q(qSEXP);
    Rcpp::traits::input_parameter< bool >::type addcom(addcomSEXP);
    rcpp_result_gen = Rcpp::wrap(FilterBySilhouetteQuantile(s, L, fallcounts, ffilcounts, falldissim, ffildissim, q, addcom));
    return rcpp_result_gen;
END_RCPP
}
// FilterBySilhouetteThreshold
Rcpp::List FilterBySilhouetteThreshold(Rcpp::NumericVector s, Rcpp::List L, std::string fallcounts, std::string ffilcounts, std::string falldissim, std::string ffildissim, float thres, bool addcom);
RcppExport SEXP _scellpam_FilterBySilhouetteThreshold(SEXP sSEXP, SEXP LSEXP, SEXP fallcountsSEXP, SEXP ffilcountsSEXP, SEXP falldissimSEXP, SEXP ffildissimSEXP, SEXP thresSEXP, SEXP addcomSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type L(LSEXP);
    Rcpp::traits::input_parameter< std::string >::type fallcounts(fallcountsSEXP);
    Rcpp::traits::input_parameter< std::string >::type ffilcounts(ffilcountsSEXP);
    Rcpp::traits::input_parameter< std::string >::type falldissim(falldissimSEXP);
    Rcpp::traits::input_parameter< std::string >::type ffildissim(ffildissimSEXP);
    Rcpp::traits::input_parameter< float >::type thres(thresSEXP);
    Rcpp::traits::input_parameter< bool >::type addcom(addcomSEXP);
    rcpp_result_gen = Rcpp::wrap(FilterBySilhouetteThreshold(s, L, fallcounts, ffilcounts, falldissim, ffildissim, thres, addcom));
    return rcpp_result_gen;
END_RCPP
}
// ClassifAsDataFrame
Rcpp::DataFrame ClassifAsDataFrame(Rcpp::List L, std::string fdist);
RcppExport SEXP _scellpam_ClassifAsDataFrame(SEXP LSEXP, SEXP fdistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type L(LSEXP);
    Rcpp::traits::input_parameter< std::string >::type fdist(fdistSEXP);
    rcpp_result_gen = Rcpp::wrap(ClassifAsDataFrame(L, fdist));
    return rcpp_result_gen;
END_RCPP
}
// ApplyPAM
Rcpp::List ApplyPAM(std::string dissim_file, int k, std::string init_method, Rcpp::Nullable<Rcpp::NumericVector> initial_med, int max_iter, int nthreads);
RcppExport SEXP _scellpam_ApplyPAM(SEXP dissim_fileSEXP, SEXP kSEXP, SEXP init_methodSEXP, SEXP initial_medSEXP, SEXP max_iterSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type dissim_file(dissim_fileSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< std::string >::type init_method(init_methodSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type initial_med(initial_medSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(ApplyPAM(dissim_file, k, init_method, initial_med, max_iter, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// GetTD
double GetTD(Rcpp::List L, std::string dissim_file);
RcppExport SEXP _scellpam_GetTD(SEXP LSEXP, SEXP dissim_fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type L(LSEXP);
    Rcpp::traits::input_parameter< std::string >::type dissim_file(dissim_fileSEXP);
    rcpp_result_gen = Rcpp::wrap(GetTD(L, dissim_file));
    return rcpp_result_gen;
END_RCPP
}
// GetJCol
Rcpp::NumericVector GetJCol(std::string fname, int ncol);
RcppExport SEXP _scellpam_GetJCol(SEXP fnameSEXP, SEXP ncolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    rcpp_result_gen = Rcpp::wrap(GetJCol(fname, ncol));
    return rcpp_result_gen;
END_RCPP
}
// GetJManyCols
Rcpp::NumericMatrix GetJManyCols(std::string fname, Rcpp::NumericVector extcols);
RcppExport SEXP _scellpam_GetJManyCols(SEXP fnameSEXP, SEXP extcolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type extcols(extcolsSEXP);
    rcpp_result_gen = Rcpp::wrap(GetJManyCols(fname, extcols));
    return rcpp_result_gen;
END_RCPP
}
// GetJColByName
Rcpp::NumericVector GetJColByName(std::string fname, std::string colname);
RcppExport SEXP _scellpam_GetJColByName(SEXP fnameSEXP, SEXP colnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type colname(colnameSEXP);
    rcpp_result_gen = Rcpp::wrap(GetJColByName(fname, colname));
    return rcpp_result_gen;
END_RCPP
}
// GetJManyColsByNames
Rcpp::NumericMatrix GetJManyColsByNames(std::string fname, Rcpp::StringVector extcolnames);
RcppExport SEXP _scellpam_GetJManyColsByNames(SEXP fnameSEXP, SEXP extcolnamesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type extcolnames(extcolnamesSEXP);
    rcpp_result_gen = Rcpp::wrap(GetJManyColsByNames(fname, extcolnames));
    return rcpp_result_gen;
END_RCPP
}
// GetSubdiag
Rcpp::NumericVector GetSubdiag(std::string fname);
RcppExport SEXP _scellpam_GetSubdiag(SEXP fnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    rcpp_result_gen = Rcpp::wrap(GetSubdiag(fname));
    return rcpp_result_gen;
END_RCPP
}
// GetJRow
Rcpp::NumericVector GetJRow(std::string fname, int nrow);
RcppExport SEXP _scellpam_GetJRow(SEXP fnameSEXP, SEXP nrowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< int >::type nrow(nrowSEXP);
    rcpp_result_gen = Rcpp::wrap(GetJRow(fname, nrow));
    return rcpp_result_gen;
END_RCPP
}
// GetJManyRows
Rcpp::NumericMatrix GetJManyRows(std::string fname, Rcpp::NumericVector extrows);
RcppExport SEXP _scellpam_GetJManyRows(SEXP fnameSEXP, SEXP extrowsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type extrows(extrowsSEXP);
    rcpp_result_gen = Rcpp::wrap(GetJManyRows(fname, extrows));
    return rcpp_result_gen;
END_RCPP
}
// GetJRowByName
Rcpp::NumericVector GetJRowByName(std::string fname, std::string rowname);
RcppExport SEXP _scellpam_GetJRowByName(SEXP fnameSEXP, SEXP rownameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type rowname(rownameSEXP);
    rcpp_result_gen = Rcpp::wrap(GetJRowByName(fname, rowname));
    return rcpp_result_gen;
END_RCPP
}
// GetJManyRowsByNames
Rcpp::NumericMatrix GetJManyRowsByNames(std::string fname, Rcpp::StringVector extrownames);
RcppExport SEXP _scellpam_GetJManyRowsByNames(SEXP fnameSEXP, SEXP extrownamesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type extrownames(extrownamesSEXP);
    rcpp_result_gen = Rcpp::wrap(GetJManyRowsByNames(fname, extrownames));
    return rcpp_result_gen;
END_RCPP
}
// JMatInfo
void JMatInfo(std::string fname, std::string fres);
RcppExport SEXP _scellpam_JMatInfo(SEXP fnameSEXP, SEXP fresSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type fres(fresSEXP);
    JMatInfo(fname, fres);
    return R_NilValue;
END_RCPP
}
// GetJRowNames
Rcpp::StringVector GetJRowNames(std::string fname);
RcppExport SEXP _scellpam_GetJRowNames(SEXP fnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    rcpp_result_gen = Rcpp::wrap(GetJRowNames(fname));
    return rcpp_result_gen;
END_RCPP
}
// GetJColNames
Rcpp::StringVector GetJColNames(std::string fname);
RcppExport SEXP _scellpam_GetJColNames(SEXP fnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    rcpp_result_gen = Rcpp::wrap(GetJColNames(fname));
    return rcpp_result_gen;
END_RCPP
}
// GetJNames
Rcpp::List GetJNames(std::string fname);
RcppExport SEXP _scellpam_GetJNames(SEXP fnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    rcpp_result_gen = Rcpp::wrap(GetJNames(fname));
    return rcpp_result_gen;
END_RCPP
}
// JWriteBin
void JWriteBin(Rcpp::NumericMatrix M, std::string fname, std::string dtype, std::string dmtype, std::string comment);
RcppExport SEXP _scellpam_JWriteBin(SEXP MSEXP, SEXP fnameSEXP, SEXP dtypeSEXP, SEXP dmtypeSEXP, SEXP commentSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type M(MSEXP);
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type dtype(dtypeSEXP);
    Rcpp::traits::input_parameter< std::string >::type dmtype(dmtypeSEXP);
    Rcpp::traits::input_parameter< std::string >::type comment(commentSEXP);
    JWriteBin(M, fname, dtype, dmtype, comment);
    return R_NilValue;
END_RCPP
}
// CalculateSilhouette
Rcpp::NumericVector CalculateSilhouette(Rcpp::NumericVector cl, std::string fdist, int nthreads);
RcppExport SEXP _scellpam_CalculateSilhouette(SEXP clSEXP, SEXP fdistSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type cl(clSEXP);
    Rcpp::traits::input_parameter< std::string >::type fdist(fdistSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(CalculateSilhouette(cl, fdist, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// NumSilToClusterSil
Rcpp::NumericMatrix NumSilToClusterSil(Rcpp::NumericVector cl, Rcpp::NumericVector s);
RcppExport SEXP _scellpam_NumSilToClusterSil(SEXP clSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type cl(clSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(NumSilToClusterSil(cl, s));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scellpam_dgCMatToJMat", (DL_FUNC) &_scellpam_dgCMatToJMat, 7},
    {"_scellpam_GetSeuratGroups", (DL_FUNC) &_scellpam_GetSeuratGroups, 1},
    {"_scellpam_SceToJMat", (DL_FUNC) &_scellpam_SceToJMat, 9},
    {"_scellpam_BuildAbundanceMatrix", (DL_FUNC) &_scellpam_BuildAbundanceMatrix, 3},
    {"_scellpam_CsvToJMat", (DL_FUNC) &_scellpam_CsvToJMat, 8},
    {"_scellpam_JMatToCsv", (DL_FUNC) &_scellpam_JMatToCsv, 4},
    {"_scellpam_ScellpamSetDebug", (DL_FUNC) &_scellpam_ScellpamSetDebug, 3},
    {"_scellpam_CalcAndWriteDissimilarityMatrix", (DL_FUNC) &_scellpam_CalcAndWriteDissimilarityMatrix, 6},
    {"_scellpam_FilterJMatByName", (DL_FUNC) &_scellpam_FilterJMatByName, 4},
    {"_scellpam_FilterBySilhouetteQuantile", (DL_FUNC) &_scellpam_FilterBySilhouetteQuantile, 8},
    {"_scellpam_FilterBySilhouetteThreshold", (DL_FUNC) &_scellpam_FilterBySilhouetteThreshold, 8},
    {"_scellpam_ClassifAsDataFrame", (DL_FUNC) &_scellpam_ClassifAsDataFrame, 2},
    {"_scellpam_ApplyPAM", (DL_FUNC) &_scellpam_ApplyPAM, 6},
    {"_scellpam_GetTD", (DL_FUNC) &_scellpam_GetTD, 2},
    {"_scellpam_GetJCol", (DL_FUNC) &_scellpam_GetJCol, 2},
    {"_scellpam_GetJManyCols", (DL_FUNC) &_scellpam_GetJManyCols, 2},
    {"_scellpam_GetJColByName", (DL_FUNC) &_scellpam_GetJColByName, 2},
    {"_scellpam_GetJManyColsByNames", (DL_FUNC) &_scellpam_GetJManyColsByNames, 2},
    {"_scellpam_GetSubdiag", (DL_FUNC) &_scellpam_GetSubdiag, 1},
    {"_scellpam_GetJRow", (DL_FUNC) &_scellpam_GetJRow, 2},
    {"_scellpam_GetJManyRows", (DL_FUNC) &_scellpam_GetJManyRows, 2},
    {"_scellpam_GetJRowByName", (DL_FUNC) &_scellpam_GetJRowByName, 2},
    {"_scellpam_GetJManyRowsByNames", (DL_FUNC) &_scellpam_GetJManyRowsByNames, 2},
    {"_scellpam_JMatInfo", (DL_FUNC) &_scellpam_JMatInfo, 2},
    {"_scellpam_GetJRowNames", (DL_FUNC) &_scellpam_GetJRowNames, 1},
    {"_scellpam_GetJColNames", (DL_FUNC) &_scellpam_GetJColNames, 1},
    {"_scellpam_GetJNames", (DL_FUNC) &_scellpam_GetJNames, 1},
    {"_scellpam_JWriteBin", (DL_FUNC) &_scellpam_JWriteBin, 5},
    {"_scellpam_CalculateSilhouette", (DL_FUNC) &_scellpam_CalculateSilhouette, 3},
    {"_scellpam_NumSilToClusterSil", (DL_FUNC) &_scellpam_NumSilToClusterSil, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_scellpam(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
