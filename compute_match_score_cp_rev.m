function matchscore = compute_match_score_cp_rev(wpdata, wloglike, testmodel, ID)

scores = OnWorkingModel(testmodel, wpdata, 1, ID);
scores = scores - wloglike; 

matchscore = max(scores);

