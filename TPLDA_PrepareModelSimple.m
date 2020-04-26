function PModel = TPLDA_PrepareModelSimple(model)

PModel.meanVec = model.meanVec;
PModel.F = model.F;
PModel.G = model.G;
PModel.Sigma = model.Sigma;

N_CONDITIONS = length(PModel.F);
N_INDIVS = length(PModel.G);
N_OBS_DIM = size(PModel.F{1}, 1);
N_HID_DIM = size(PModel.F{1}, 2);
N_HID_DIM_NOISE = size(PModel.G{1}{1}, 2);

for (cCond = 1:N_CONDITIONS)
    PModel.invSigmaEst{cCond} = 1./PModel.Sigma{cCond};
    PModel.SF{cCond} = PModel.F{cCond}.*repmat(PModel.invSigmaEst{cCond}, 1, N_HID_DIM);
    PModel.LogDetSigma{cCond} = sum(log(PModel.Sigma{cCond}));

    for (cInd = 1:N_INDIVS)
        PModel.SG{cInd}{cCond} = PModel.G{cInd}{cCond}.*repmat(PModel.invSigmaEst{cCond}, 1, N_HID_DIM_NOISE);
        PModel.A{cInd}{cCond} = [PModel.F{cCond} PModel.G{cInd}{cCond}];
        PModel.SA{cInd}{cCond} = [PModel.SF{cCond} PModel.SG{cInd}{cCond}];

        PModel.IATSA{cInd}{cCond} = PModel.A{cInd}{cCond}'*PModel.SA{cInd}{cCond} + eye(N_HID_DIM + N_HID_DIM_NOISE, N_HID_DIM + N_HID_DIM_NOISE);
        PModel.InvIATSA{cInd}{cCond} = inv(PModel.IATSA{cInd}{cCond});
        PModel.LogDetIATSA{cInd}{cCond} = logdet(PModel.IATSA{cInd}{cCond}, 'chol');

        PModel.FrLogTerm{cInd}{cCond} = -N_OBS_DIM/2*log(2*pi) - 0.5*PModel.LogDetSigma{cCond} - 0.5*PModel.LogDetIATSA{cInd}{cCond}; 
    end;
end;

