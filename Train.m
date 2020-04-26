function trainresults = Train()

FEATURE_FOLDER = '../Features1A';
load([FEATURE_FOLDER '/lbp_faces_ytc_3.mat']);
load([FEATURE_FOLDER '/splits_ytc.mat'])

NUM_SUBJECTS = numel(lbp_faces);
NUM_SAMPLES = 3;

F_FACTORS = [42];
G_FACTORS = [42];
N_ITER = 16;

for ctrial = 1:10
    [allData allIDs allConditions] = LoadData(lbp_faces, splits, NUM_SUBJECTS, NUM_SAMPLES, ctrial);

    for cfac = 1:length(F_FACTORS)
        model = TPLDA_Train_Hetero(allData, allIDs, allConditions, F_FACTORS(cfac), G_FACTORS(cfac), N_ITER);

        modelname = sprintf('%d_model_%d_%d.mat', ctrial, F_FACTORS(cfac), G_FACTORS(cfac));
        save(modelname, 'model');
    end;
end;

%===============================================================================

function [allData allIDs allConditions] = LoadData(fets, spts, numsubjects, numsamples, trial)

allData = [];
allIDs = [];
allConditions = [];

cp = 0;
cfold = spts{trial}.train;

for i = 1:numsubjects
    fprintf('load subject-%d\n', i);
        
    for j = 1:numsamples
        cp = cp + 1;
		cinds = cfold{i}{j};
        
        cfeat = im2double(fets{i}{cinds(1)}{cinds(2)});
        nfeat = size(cfeat, 2);

        allData = [allData cfeat];

        cnum = size(allIDs, 1);
        trmat = zeros(cnum, 1);
        
        bmat = zeros(nfeat, cp);
        bmat(:, cp) = 1;
        allIDs = [allIDs trmat; bmat];

        cconds = ones(nfeat, 1);
        allConditions = [allConditions; cconds];
    end;
end;

