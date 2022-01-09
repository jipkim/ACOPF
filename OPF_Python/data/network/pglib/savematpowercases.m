testsystemlist = ... 
   ["pglib_opf_case3_lmbd",
    "pglib_opf_case5_pjm",
    "pglib_opf_case14_ieee",
    "pglib_opf_case24_ieee_rts",
    "pglib_opf_case30_as",
    "pglib_opf_case30_ieee",
    "pglib_opf_case39_epri",
    "pglib_opf_case57_ieee",
    "pglib_opf_case57_ieee",
    "pglib_opf_case73_ieee_rts",
    "pglib_opf_case89_pegase",
    "pglib_opf_case118_ieee",
    "pglib_opf_case162_ieee_dtc",
    "pglib_opf_case179_goc",
    "pglib_opf_case200_activ",
    "pglib_opf_case240_pserc",
    "pglib_opf_case300_ieee",
    "pglib_opf_case500_goc",
    "pglib_opf_case588_sdet",
    "pglib_opf_case793_goc",
    "pglib_opf_case1354_pegase",
    "pglib_opf_case1888_rte"]


for ii = 1:length(testsystemlist)
    clear mpc;
    testsystem = '/Users/jipkim/Dropbox/Matlab/pglib-opf/' + testsystemlist(ii) + '.m';
    filename = testsystemlist(ii) + '.mat';
    mpc = loadcase(testsystem);
    save(filename, 'mpc') ;
end
%%
testsystem = 'case9';
filename = [testsystem, '.mat'];
mpc = loadcase(testsystem);
save(filename, mpc) ;
