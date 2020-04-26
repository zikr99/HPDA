function matchscore = compute_match_score(testdata, wmodel, prepmodel, cind)

N_DATA = size(testdata, 2);

wpdata = TPLDA_PrepareData(testdata, ones(N_DATA, 1), prepmodel, cind);
mmodel = IncWorkingModel(wmodel, wpdata, ones(N_DATA, 1), cind);

matchscore = mmodel.Loglike - wmodel.Loglike;

