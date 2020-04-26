function rmodel = TPLDA_Train_Hetero(data, imageIDs, conditions, N_FAC, N_G, N_ITER)

%define number of factors
fprintf('Estimating %d factors, %d g factors, %d iterations\n', N_FAC, N_G, N_ITER);	

OBS_DIM = size(data, 1);
N_DATA = size(data, 2);
N_INDIVS = size(imageIDs, 2);

%condition is a vector that states which of several conditions each data point is in.
nCondition = max(conditions);

allclusters = [];

for (cCond = 1:nCondition)
    %calculate means and remove from the dataset
    condIndices = find(conditions == cCond);

    %meanVec{cCond} = mean(data(:, condIndices), 2);
    meanVec{cCond} = zeros(OBS_DIM, 1);

    for i = condIndices'
        meanVec{cCond} = meanVec{cCond} + data(:, i);
    end;
    
    meanVec{cCond} = meanVec{cCond}/length(condIndices);

    %data(:, condIndices) = data(:, condIndices) - repmat(meanVec{cCond}, 1, length(condIndices));
    for i = condIndices'
        data(:, i) = data(:, i) - meanVec{cCond};
    end;

    %clusters{cCond} = (data(:, condIndices)*imageIDs(condIndices, :))*(diag(1./(sum(imageIDs(condIndices, :))))); 
    clusters = zeros(OBS_DIM, N_INDIVS);

    for i = 1:N_INDIVS
        for j = condIndices'
            clusters(:, i) = clusters(:, i) + data(:, j)*imageIDs(j, i); 
        end;
 
        clusters(:, i) = clusters(:, i)/sum(imageIDs(condIndices, i));        
    end;

    intraIndices = [];

    for (cCi = condIndices')
        intraIdx = find(imageIDs(cCi, :));
        intraIndices = [intraIndices intraIdx];
    end;

    %intras{cCond} = data(:, condIndices) - clusters{cCond}(:, intraIndices);
    intras = zeros(OBS_DIM, length(condIndices)); 

    for i = 1:length(condIndices)
        intras(:, i) = data(:, condIndices(i)) - clusters(:, intraIndices(i));
    end;

    %[inevecs inssq{cCond}] = PCA_Train2(intras, N_G);

    for i = 1:N_INDIVS
        %GEst{i}{cCond} = inevecs; 
        GEst{i}{cCond} = randn(OBS_DIM, N_G);
    end;
    
    allclusters = [allclusters; clusters];
end; 

%[evecs ssq] = PCA_Train2(allclusters, N_FAC); 

%initialize factors and noises
for (cCond = 1:nCondition)
    %FEst{cCond} = evecs((cCond - 1)*OBS_DIM + 1:cCond*OBS_DIM, :);
    %SigmaEst{cCond} = (ssq + inssq{cCond})*ones(OBS_DIM, 1);
    FEst{cCond} = randn(OBS_DIM, N_FAC);
    SigmaEst{cCond} = 0.01*ones(OBS_DIM, 1);
end;

%model.meanVec = meanVec;
%model.F = FEst;
%model.G = GEst;
%model.Sigma = SigmaEst;

%modelname = sprintf('model_%d_%d_%d.mat', N_FAC, N_G, 0);
%save(modelname, 'model');

%for each iteration
for (cIter = 1:N_ITER)
    fprintf('Iter %d\n', cIter);

    %calculate expected value of h and hth
    [Eh EhhSum EhhSumInd] = getExpectedValuesTied(FEst, GEst, SigmaEst, data, imageIDs, conditions);

    for (cCond = 1:nCondition)
        condIndices = find(conditions == cCond);

        xhSum = zeros(OBS_DIM, N_FAC);
        
        for (cData = condIndices')
            xhSum = xhSum + data(:, cData)*Eh(1:N_FAC, cData)';  
        end;

        scSum = zeros(OBS_DIM, N_FAC);
  
        for (i = 1:N_INDIVS)
            scSum = scSum + GEst{i}{cCond}*EhhSumInd{i}{cCond}(N_FAC + 1:N_FAC + N_G, 1:N_FAC);
        end;
    
        %update F
        FEst{cCond} = (xhSum - scSum)*inv(EhhSum{cCond}(1:N_FAC, 1:N_FAC));

        for (i = 1:N_INDIVS)
            xhSumInd{i} = zeros(OBS_DIM, N_G);
        end;

        for (cData = condIndices')
            inInd = find(imageIDs(cData, :));  

            xhSumInd{inInd} = xhSumInd{inInd} + data(:, cData)*Eh(N_FAC + 1:N_FAC + N_G, cData)';  
        end;

        for (i = 1:N_INDIVS)
            GEst{i}{cCond} = (xhSumInd{i} - FEst{cCond}*EhhSumInd{i}{cCond}(1:N_FAC, N_FAC + 1:N_FAC + N_G))*inv(EhhSumInd{i}{cCond}(N_FAC + 1:N_FAC + N_G, N_FAC + 1:N_FAC + N_G));    
        end;

        %update Sigma 
        %SigmaEst{cCond} = mean(data(:, condIndices).*data(:, condIndices) - (FG{cCond}*Eh(:, condIndices)).*data(:, condIndices), 2); 
        SigmaEst{cCond} = zeros(OBS_DIM, 1);  

        for (i = condIndices')
            inInd = find(imageIDs(i, :)); 

            SigmaEst{cCond} = SigmaEst{cCond} + data(:, i).*data(:, i) - ([FEst{cCond} GEst{inInd}{cCond}]*Eh(:, i)).*data(:, i);
        end;

        SigmaEst{cCond} = SigmaEst{cCond}/length(condIndices);
    end;

    %model.meanVec = meanVec;
    %model.F = FEst;
    %model.G = GEst;
    %model.Sigma = SigmaEst;

    %modelname = sprintf('model_%d_%d_%d.mat', N_FAC, N_G, cIter);
    %save(modelname, 'model');
end;

%for (cCond = 1:nCondition)
    %condIndex = find(condition == cCond);

    %deal with scale issues - deals with slow convergence to final answer.
    %EhSD = sqrt(diag(cov(Eh(:, condIndex)')));
    %FG{cCond} = FG{cCond}.*repmat(EhSD', OBS_DIM, 1);
     
    %FEst{cCond} = FG{cCond}(:, 1:N_FAC_EST);
    %GEst{cCond} = FG{cCond}(:, N_FAC_EST + 1:end);

    %sigmaEst{cCond} = sigmaEst{cCond}*100;
%end;

rmodel.meanVec = meanVec;
rmodel.F = FEst;
rmodel.G = GEst;
rmodel.Sigma = SigmaEst;

%=====================================================================================

function [Eh, EhhSum, EhhSumInd] = getExpectedValuesTied(FEst, GEst, SigmaEst, x, imageIDs, conditions)

N_DATA = size(x, 2);
N_INDIVS = size(imageIDs, 2);
N_CONDITIONS = max(conditions);

N_HID_DIM = size(FEst{1}, 2);
N_HID_DIM_NOISE = size(GEst{1}{1}, 2);

for (cCond = 1:N_CONDITIONS)
    invSigmaEst{cCond} = 1./SigmaEst{cCond};
    FTSF{cCond} = FEst{cCond}'*(FEst{cCond}.*repmat(invSigmaEst{cCond}, 1, N_HID_DIM));    

    for (i = 1:N_INDIVS)   
        FTSG{i}{cCond} = FEst{cCond}'*(GEst{i}{cCond}.*repmat(invSigmaEst{cCond}, 1, N_HID_DIM_NOISE));
        GTSF{i}{cCond} = GEst{i}{cCond}'*(FEst{cCond}.*repmat(invSigmaEst{cCond}, 1, N_HID_DIM));
        GTSG{i}{cCond} = GEst{i}{cCond}'*(GEst{i}{cCond}.*repmat(invSigmaEst{cCond}, 1, N_HID_DIM_NOISE));  

        GTSGPlusI{i}{cCond} = GTSG{i}{cCond} + eye(N_HID_DIM_NOISE, N_HID_DIM_NOISE);
	  InvGTSGPlusI{i}{cCond} = inv(GTSGPlusI{i}{cCond});
    	  DetGTSGPlusI{i}{cCond} = det(GTSGPlusI{i}{cCond}); 

        BDI{i}{cCond} = FTSG{i}{cCond}*InvGTSGPlusI{i}{cCond};
        DIC{i}{cCond} = InvGTSGPlusI{i}{cCond}*GTSF{i}{cCond};
        BDICi{i}{cCond} = BDI{i}{cCond}*GTSF{i}{cCond};
        AiMinBDICi{i}{cCond} = FTSF{cCond} - BDICi{i}{cCond};      
    end;
end;

Eh = zeros(N_HID_DIM + N_HID_DIM_NOISE, N_DATA);

for (cCond = 1:N_CONDITIONS)
    EhhSum{cCond}  = zeros(N_HID_DIM + N_HID_DIM_NOISE, N_HID_DIM + N_HID_DIM_NOISE);

    for (i = 1:N_INDIVS)
        EhhSumInd{i}{cCond} = zeros(N_HID_DIM + N_HID_DIM_NOISE, N_HID_DIM + N_HID_DIM_NOISE);
    end;
end;
 
for (cInd = 1:N_INDIVS)
    %figure out how many data points we are combining here
    imIndices = find(imageIDs(:, cInd));
    nImg = length(imIndices);

    X00 = zeros(N_HID_DIM, N_HID_DIM);

    for (i = 1:nImg)
        X00 = X00 + AiMinBDICi{cInd}{conditions(imIndices(i))};
    end;

    X11 = inv(X00 + eye(N_HID_DIM, N_HID_DIM));

    for (cCond = 1:N_CONDITIONS)
        X12{cCond} = -X11*BDI{cInd}{cCond};
        X21{cCond} = -DIC{cInd}{cCond}*X11;
    end;

    for (cCond1 = 1:N_CONDITIONS)
        for (cCond2 = 1:N_CONDITIONS)
            X22mDI{cCond1}{cCond2} = -X21{cCond1}*BDI{cInd}{cCond2};
        end;
    end;

    SFTSX = zeros(N_HID_DIM, 1);

    for (cIm = 1:length(imIndices))
        SFTSX = SFTSX + FEst{conditions(imIndices(cIm))}'*(x(:, imIndices(cIm)).*invSigmaEst{conditions(imIndices(cIm))});
    end;

    for (cIm = 1:length(imIndices))
        GTSX{cIm} = GEst{cInd}{conditions(imIndices(cIm))}'*(x(:, imIndices(cIm)).*invSigmaEst{conditions(imIndices(cIm))});
    end;

    Ehi = X11*SFTSX;

    for (cIm = 1:length(imIndices))
        Ehi = Ehi + X12{conditions(imIndices(cIm))}*GTSX{cIm};
    end;

    for (cIm1 = 1:length(imIndices))
        Ewi = X21{conditions(imIndices(cIm1))}*SFTSX;

        for (cIm2 = 1:length(imIndices))
            Ewi = Ewi + X22mDI{conditions(imIndices(cIm1))}{conditions(imIndices(cIm2))}*GTSX{cIm2};
        end;

        Ewi = Ewi + InvGTSGPlusI{cInd}{conditions(imIndices(cIm1))}*GTSX{cIm1};

        Eh(:, imIndices(cIm1)) = [Ehi; Ewi];

        EhhSum{conditions(imIndices(cIm1))} = EhhSum{conditions(imIndices(cIm1))} + ... 
            [X11 X12{conditions(imIndices(cIm1))}; X21{conditions(imIndices(cIm1))} X22mDI{conditions(imIndices(cIm1))}{conditions(imIndices(cIm1))} + ...
            InvGTSGPlusI{cInd}{conditions(imIndices(cIm1))}] + [Ehi; Ewi]*[Ehi; Ewi]';

        EhhSumInd{cInd}{conditions(imIndices(cIm1))} = EhhSumInd{cInd}{conditions(imIndices(cIm1))} + ...
            [X11 X12{conditions(imIndices(cIm1))}; X21{conditions(imIndices(cIm1))} X22mDI{conditions(imIndices(cIm1))}{conditions(imIndices(cIm1))} + ...
            InvGTSGPlusI{cInd}{conditions(imIndices(cIm1))}] + [Ehi; Ewi]*[Ehi; Ewi]';
    end;
end;

%==========================================================================

function drawPictures(pixelData,FEst,sigmaEst,meanVec,imageID)

%display standard deviation image
sigmaIm = sqrt(sigmaEst);
sigmaIm = reshape(sigmaIm,sqrt(length(sigmaIm)),sqrt(length(sigmaIm)));
figure; set(gcf,'Color',[1 1 1]); imagesc(sigmaIm); axis off; axis image; colormap(gray);
set(gcf,'Name','Standard Deviation Image');

%display first factors
for (cComp = 1:5)
    figure; set(gcf,'Color',[1 1 1]); 
    factor = FEst(:,cComp);
    factor = reshape(factor,sqrt(length(factor)),sqrt(length(factor)));
    imagesc(factor); axis off; axis image; colormap(gray);
    set(gcf,'Name',sprintf('Factor %d',cComp));
end;

observed = pixelData-repmat(meanVec,1,size(pixelData,2));
[Eh Ehh] = getExpectedValues(FEst,sigmaEst,observed,imageID);

%show reconstruction so we can check we haven't screwed up...
for (cIm = 1:4:20)
    figure;set(gcf,'Color',[1 1 1]);set(gcf,'Position',[358   433   560   198]);
    set(gcf,'Name',sprintf('Reconstructed Image %d and Original Image %d',cIm,cIm));

    recon = FEst*Eh(:,cIm);
    recon = recon+meanVec;   
    recon = reshape(recon,sqrt(length(recon)),sqrt(length(recon)));
    subplot(1,2,1); imagesc(recon,[0 255]);colormap(gray); axis off; axis image; colorbar;
    im1Orig = observed(:,cIm)+meanVec;
    orig = reshape(im1Orig,sqrt(length(im1Orig)),sqrt(length(im1Orig)));
    subplot(1,2,2); imagesc(orig,[0 255]);colormap(gray); axis off; axis image; colorbar;
end;

