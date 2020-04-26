function [frscores frids vdscores vdids] = Test(trial)

MODEL_FOLDER = '../Models1';
FEATURE_FOLDER = '../Features1';

F_FACTORS = [32]; %[1 2 4 8 16 24 32 42];
G_FACTORS = [32]; %[1 2 4 8 16 24 32 42];

allDataTrain = [];
allIDsTrain = [];
allDataTest = [];
allIDsTest = [];

[allDataTrain allIDsTrain] = LoadDataTrain(trial);
[allDataTest allIDsTest] = LoadDataTest(trial);

%featfile = sprintf('%s/Test%d.mat', FEATURE_FOLDER, trial);
%fprintf('load %s\n', featfile);
%load(featfile);

testresults = [];

for cfac = 1:length(F_FACTORS)
    for cIter = 0:0
    modelname = sprintf('%s/%d_model_%d_%d.mat', MODEL_FOLDER, trial, F_FACTORS(cfac), G_FACTORS(cfac));
    %modelname = sprintf('%s/model_%d_%d_%d.mat', MODEL_FOLDER, F_FACTORS(cfac), G_FACTORS(cfac), cIter);
    fprintf('Load model: %s\n', modelname);

    load(modelname);
    TPLDAModel = model;

    prepmodel = TPLDA_PrepareModel(TPLDAModel);
    simpmodel = TPLDA_PrepareModelSimple(TPLDAModel);
    simpmodel12 = TPLDA_PrepareModelSimple12(TPLDAModel);

    fprintf('Preparation: \n');

    wpdata = [];
    wmodels = [];
    wlls = [];

    numindivs = length(allDataTrain);

    for cind = 1:numindivs
        fprintf('train-%d ', cind);

        ndata = size(allDataTrain{cind}, 2);
        ndata4 = floor(ndata/4) + 1;
        ndata2 = floor(ndata/2) + 1;

        mndata = mean(allDataTrain{cind}, 2);
        
        dd = allDataTrain{cind} - repmat(mndata, 1, ndata);
        dm = dd.*dd;        
        dn = sum(dm, 1);

        [se oi] = sort(dn);

        wpdata{cind} = TPLDA_PrepareData(allDataTrain{cind}(:, oi(1:ndata4)), ones(ndata4, 1), prepmodel, cind);
        wmodels{cind} = InitWorkingModel(prepmodel, wpdata{cind}, ones(ndata4, 1), cind);
        wlls{cind} = OnSimpleModel(simpmodel, wpdata{cind}, 1, cind);

        wpdatamn{cind} = TPLDA_PrepareData(allDataTrain{cind}(:, oi(1:ndata2)), ones(ndata2, 1), prepmodel, cind);
        wmodelsmn{cind} = InitWorkingModel(prepmodel, wpdatamn{cind}, ones(ndata2, 1), cind);
        wllsmn{cind} = OnSimpleModel(simpmodel, wpdatamn{cind}, 1, cind);
    end;

    fprintf('\n');
    fprintf('Identification: \n');

    numsubjects = length(allDataTest);
    %numcorrects = 0;
    %numall = 0;
    %indres = [];

    matscores1 = [];
    matscores2 = [];
    matscores3 = [];
    matids = [];

    matscs1 = [];
    matscs2 = [];
    matscs3 = [];
    matds = [];     

    for csub = 1:numsubjects
        fprintf('%d ', csub);   
        %numcorrects = 0;
        %numall = 0;

        for cfr = 1:size(allDataTest{csub}, 2)     
            cs1 = [];            
            cs2 = [];
            cs3 = []; 

            for cind = 1:numindivs
                cscore1 = compute_match_score_cp(allDataTest{csub}(:, cfr), wpdata{cind}, wlls{cind}, prepmodel, simpmodel12, cind);
                cscore2 = compute_match_score_cp(allDataTest{csub}(:, cfr), wpdatamn{cind}, wllsmn{cind}, prepmodel, simpmodel12, cind);
                cscore3 = compute_match_score(allDataTest{csub}(:, cfr), wmodels{cind}, prepmodel, cind);

                cs1 = [cs1; cscore1];
                cs2 = [cs2; cscore2];
                cs3 = [cs3; cscore3];
            end;

            %[mv pos] = max(matscores);

            %if (allIDsTest{csub} == allIDsTrain{pos})
            %    numcorrects = numcorrects + 1;
            %end; 

            %numall = numall + 1;

            matscores1 = [matscores1 cs1];            
            matscores2 = [matscores2 cs2];
            matscores3 = [matscores3 cs3];
            matids = [matids csub];
        end;

        cs1 = [];            
        cs2 = [];
        cs3 = []; 

        for cind = 1:numindivs
            ntestdata = size(allDataTest{csub}, 2);
            wptestdata = TPLDA_PrepareData(allDataTest{csub}, ones(ntestdata, 1), prepmodel, cind);
            wtestmodel = InitWorkingModel(prepmodel, wptestdata, ones(ntestdata, 1), cind);
            
            cscore1 = compute_match_score_cp_rev(wpdata{cind}, wlls{cind}, wtestmodel, cind);
            cscore2 = compute_match_score(allDataTest{csub}, wmodelsmn{cind}, prepmodel, cind);
            cscore3 = compute_match_score(allDataTest{csub}, wmodels{cind}, prepmodel, cind);

            cs1 = [cs1; cscore1];
            cs2 = [cs2; cscore2];
            cs3 = [cs3; cscore3];
        end;

        matscs1 = [matscs1 cs1];            
        matscs2 = [matscs2 cs2];
        matscs3 = [matscs3 cs3];
        matds = [matds csub];

        %indres = [indres; numcorrects/numall];
    end;

    fprintf('\n');

    %disp(numcorrects/numall);
    %testresults = [testresults; numcorrects/numall];
    %testresults{cfac} = indres;

    frscores{cfac}{1} = matscores1;
    frscores{cfac}{2} = matscores2;
    frscores{cfac}{3} = matscores3;
    frids{cfac} = matids;
    
    vdscores{cfac}{1} = matscs1;
    vdscores{cfac}{2} = matscs2;
    vdscores{cfac}{3} = matscs3;
    vdids{cfac} = matds;
    end;
end;

%==============================================================================

function [allData allIDs] = LoadDataTrain(trial)

datafolder = '../Features1';
numsubjects = 43;

allData = [];
allIDs = [];
allConditions = [];

for i = 0:(numsubjects - 1)
    featfile = sprintf('%s/%d_train%d.txt', datafolder, i, trial);
    fprintf('load %s\n', featfile);

    cfeat = load(featfile);
    nfeat = size(cfeat, 1);

    allData{i + 1} = cfeat';
    allIDs{i + 1} = i + 1;
end;

%==============================================================================

function [allData allIDs] = LoadDataTest(trial)

datafolder = '../Features1';
numsubjects = 43;

allData = [];
allIDs = [];
allConditions = [];

for i = 0:(numsubjects - 1)
    featfile = sprintf('%s/%d_test%d.txt', datafolder, i, trial);
    fprintf('load %s\n', featfile);

    cfeat = load(featfile);
    nfeat = size(cfeat, 1);

    allData{i + 1} = cfeat';
    allIDs{i + 1} = i + 1;
end;

