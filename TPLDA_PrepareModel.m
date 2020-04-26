function PModel = TPLDA_PrepareModel(model)

PModel.meanVec = model.meanVec;
PModel.F = model.F;
PModel.G = model.G;
PModel.Sigma = model.Sigma;

N_CONDITIONS = length(PModel.F);
N_INDIVS = length(PModel.G);
N_HID_DIM = size(PModel.F{1}, 2);
N_HID_DIM_NOISE = size(PModel.G{1}{1}, 2);

for (cCond = 1:N_CONDITIONS)
    PModel.invSigmaEst{cCond} = 1./PModel.Sigma{cCond};
    PModel.FTSF{cCond} = PModel.F{cCond}'*(PModel.F{cCond}.*repmat(PModel.invSigmaEst{cCond}, 1, N_HID_DIM));

    for (cInd = 1:N_INDIVS)
        PModel.FTSG{cInd}{cCond} = PModel.F{cCond}'*(PModel.G{cInd}{cCond}.*repmat(PModel.invSigmaEst{cCond}, 1, N_HID_DIM_NOISE));
        PModel.GTSF{cInd}{cCond} = PModel.G{cInd}{cCond}'*(PModel.F{cCond}.*repmat(PModel.invSigmaEst{cCond}, 1, N_HID_DIM));
        PModel.GTSG{cInd}{cCond} = PModel.G{cInd}{cCond}'*(PModel.G{cInd}{cCond}.*repmat(PModel.invSigmaEst{cCond}, 1, N_HID_DIM_NOISE));  

        PModel.GTSGPlusI{cInd}{cCond} = PModel.GTSG{cInd}{cCond} + eye(N_HID_DIM_NOISE, N_HID_DIM_NOISE);
        PModel.InvGTSGPlusI{cInd}{cCond} = inv(PModel.GTSGPlusI{cInd}{cCond});
        PModel.LogDetGTSGPlusI{cInd}{cCond} = logdet(PModel.GTSGPlusI{cInd}{cCond}, 'chol'); 

        PModel.BDI{cInd}{cCond} = PModel.FTSG{cInd}{cCond}*PModel.InvGTSGPlusI{cInd}{cCond};
        PModel.DIC{cInd}{cCond} = PModel.InvGTSGPlusI{cInd}{cCond}*PModel.GTSF{cInd}{cCond};
        PModel.BDICi{cInd}{cCond} = PModel.BDI{cInd}{cCond}*PModel.GTSF{cInd}{cCond};
        PModel.AiMinBDICi{cInd}{cCond} = PModel.FTSF{cCond} - PModel.BDICi{cInd}{cCond};   
    end;
end;

