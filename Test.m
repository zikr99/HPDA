function [fs fi vs vi galids] = Test(trial)

MODEL_FOLDER = '../Models1A';
FEATURE_FOLDER = '../Features1A';

F_FACTORS = [42];
G_FACTORS = [42];

allDataTrain = [];
allIDsTrain = [];
allDataTest = [];
allIDsTest = [];

load([FEATURE_FOLDER '/lbp_faces_ytc_3.mat']);
load([FEATURE_FOLDER '/splits_ytc.mat']);

numsubs = numel(lbp_faces);
numtrain = 3;
numtest = 6;

[allDataTrain allIDsTrain] = LoadDataTrain(lbp_faces, splits, trial, numsubs, numtrain);
[allDataTest allIDsTest] = LoadDataTest(lbp_faces, splits, trial, numsubs, numtest);

for cfac = 1:length(F_FACTORS)
    modelname = sprintf('%s/%d_model_%d_%d.mat', MODEL_FOLDER, trial, F_FACTORS(cfac), G_FACTORS(cfac));
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
        wpdata{cind} = TPLDA_PrepareData(allDataTrain{cind}, ones(ndata, 1), prepmodel, cind);
        wmodels{cind} = InitWorkingModel(prepmodel, wpdata{cind}, ones(ndata, 1), cind);
        wlls{cind} = OnSimpleModel(simpmodel, wpdata{cind}, 1, cind);

        mndata = mean(allDataTrain{cind}, 2);
        wpdatamn{cind} = TPLDA_PrepareData(mndata, [1], prepmodel, cind);
        wmodelsmn{cind} = InitWorkingModel(prepmodel, wpdatamn{cind}, [1], cind);
        wllsmn{cind} = OnSimpleModel(simpmodel, wpdatamn{cind}, 1, cind);
    end;

    fprintf('\n');
    fprintf('Identification: \n');

    numsubjects = length(allDataTest);

    matscores1 = [];
    matscores2 = [];
    matids = [];

    matscs1 = [];
    matscs2 = [];
    matscs3 = [];
    matscs4 = [];
    matscs5 = [];
    matds = [];     

    for csub = 1:numsubjects
        fprintf('%d ', csub); 

		ntestdata = size(allDataTest{csub}, 2);

        cs1 = [];            
        cs2 = [];

        for cind = 1:numindivs
            cscore1 = compute_match_score_cp(allDataTest{csub}, wpdata{cind}, wlls{cind}, prepmodel, simpmodel12, cind);
            cscore2 = compute_match_score_cp(allDataTest{csub}, wpdatamn{cind}, wllsmn{cind}, prepmodel, simpmodel12, cind);
            
            cs1 = [cs1; cscore1];
            cs2 = [cs2; cscore2];
        end;

        matscores1 = [matscores1 cs1];            
        matscores2 = [matscores2 cs2];
        matids = [matids allIDsTest{csub}*ones(1, ntestdata)];

        cs1 = sum(cs1, 2);
        cs2 = sum(cs2, 2);
        cs3 = [];            
        cs4 = [];
        cs5 = []; 

        for cind = 1:numindivs
            wptestdata = TPLDA_PrepareData(allDataTest{csub}, ones(ntestdata, 1), prepmodel, cind);
            wtestmodel = InitWorkingModel(prepmodel, wptestdata, ones(ntestdata, 1), cind);
            
            cscore3 = compute_match_score(allDataTest{csub}, wmodels{cind}, prepmodel, cind);
            cscore4 = compute_match_score(allDataTest{csub}, wmodelsmn{cind}, prepmodel, cind);
            cscore5 = compute_match_score_cp_rev(wpdata{cind}, wlls{cind}, wtestmodel, cind);

            cs3 = [cs3; cscore3];
            cs4 = [cs4; cscore4];
            cs5 = [cs5; cscore5];
        end;

        matscs1 = [matscs1 cs1];            
        matscs2 = [matscs2 cs2];
        matscs3 = [matscs3 cs3];
        matscs4 = [matscs4 cs4];
        matscs5 = [matscs5 cs5];
        matds = [matds allIDsTest{csub}];
    end;

    fprintf('\n');
    
    fs{cfac}{1} = matscores1;
    fs{cfac}{2} = matscores2;
    fi{cfac} = matids;
    
    vs{cfac}{1} = matscs1;
    vs{cfac}{2} = matscs2;
    vs{cfac}{3} = matscs3;
    vs{cfac}{4} = matscs4;
    vs{cfac}{5} = matscs5;
    vi{cfac} = matds;
    
    galids{cfac} = [];
    
    for cind = 1:numindivs
        galids{cfac} = [galids{cfac}; allIDsTrain{cind}];
    end;
    
    outfile = sprintf('%d_results_%d_%d.mat', trial, F_FACTORS(cfac), G_FACTORS(cfac));
    save(outfile, 'fs', 'fi', 'vs', 'vi', 'galids');
end;

%==============================================================================

function [allData allIDs] = LoadDataTrain(fets, spts, trial, numsubjects, numsamples)

allData = [];
allIDs = [];
allConditions = [];

cp = 0;
cfold = spts{trial}.train;

for i = 1:numsubjects
    fprintf('load subject-train-%d\n', i);
    
    for j = 1:numsamples
        cp = cp + 1;
		cinds = cfold{i}{j};
        
        cfeat = im2double(fets{i}{cinds(1)}{cinds(2)});

        allData{cp} = cfeat;
        allIDs{cp} = i;
    end;
end;

%==============================================================================

function [allData allIDs] = LoadDataTest(fets, spts, trial, numsubjects, numsamples)

allData = [];
allIDs = [];
allConditions = [];

cp = 0;
cfold = spts{trial}.test;

for i = 1:numsubjects
    fprintf('load subject-test-%d\n', i);
    
    for j = 1:numsamples
        cp = cp + 1;
		cinds = cfold{i}{j};
    
        cfeat = im2double(fets{i}{cinds(1)}{cinds(2)});
        nfeat = size(cfeat, 2);
        
        mstep = floor(nfeat/10);
        
        if mstep == 0
            mstep = 1;
        end;
        
        cfeat = cfeat(:, 1:mstep:nfeat);

        allData{cp} = cfeat;
        allIDs{cp} = i;
    end;
end;


