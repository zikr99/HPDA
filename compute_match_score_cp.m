function matchscore = compute_match_score_cp(testdata, wpdata, wloglike, ppldamodel, spmodel, ID)

N_DATA1 = size(testdata, 2);
N_DATA2 = size(wpdata.XM, 2);

wptestdata = TPLDA_PrepareData(testdata, ones(N_DATA1, 1), ppldamodel, ID);
scores = OnSimpleModel12(spmodel, wptestdata, wpdata, 1, ID);

for c1 = 1:N_DATA1
    scores(:, c1) = scores(:, c1) - wloglike';
end;

if (N_DATA2 > 1) 
    matchscore = max(scores);
else
    matchscore = scores;
end;

