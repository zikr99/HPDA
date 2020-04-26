function NModel = IncWorkingModel(WModel, PData, PConditions, ID)

N_CONDITIONS = length(WModel.F);
N_HID_DIM = size(WModel.F{1}, 2);
N_HID_DIM_NOISE = size(WModel.G{1}{1}, 2);
N_OBS_DIM = size(PData.XM, 1);
N_DATA = size(PData.XM, 2);

NModel = WModel;

for (i = 1:N_DATA)
    NModel.X00 = NModel.X00 + NModel.AiMinBDICi{ID}{PConditions(i)};
end;

NModel.X11 = inv(NModel.X00 + eye(N_HID_DIM, N_HID_DIM));

for (cCond = 1:N_CONDITIONS)
    NModel.X12{cCond} = -NModel.X11*NModel.BDI{ID}{cCond};
    NModel.X21{cCond} = -NModel.DIC{ID}{cCond}*NModel.X11;
end;

for (cCond1 = 1:N_CONDITIONS)
    for (cCond2 = 1:N_CONDITIONS)
        NModel.X22mDI{cCond1}{cCond2} = -NModel.X21{cCond1}*NModel.BDI{ID}{cCond2};
    end;
end;

NModel.LogScDet1 = logdet(NModel.X00);

for (i = 1:N_DATA)
    NModel.LogFrDet = NModel.LogFrDet + sum(log(NModel.Sigma{PConditions(i)}));
    NModel.LogScDet2 = NModel.LogScDet2 + NModel.LogDetGTSGPlusI{ID}{PConditions(i)};
end;

for (i = 1:N_DATA)
    NModel.SFTSX = NModel.SFTSX + PData.FTSX(:, i);
end;

for (cCond = 1:N_CONDITIONS)
    iIndices = find(PConditions == cCond);

    for (i = iIndices')
        NModel.SGTSX{cCond} = NModel.SGTSX{cCond} + PData.GTSX(:, i);  
    end;
end;

for (i = 1:N_DATA)
    NModel.SXTSX = NModel.SXTSX + PData.XTSX(i);
    NModel.SGXGI = NModel.SGXGI + PData.GXGI(i);
end;

NModel.V11 = NModel.SFTSX'*NModel.X11*NModel.SFTSX;

for (cCond = 1:N_CONDITIONS)
    NModel.V12{cCond} = NModel.SFTSX'*NModel.X12{cCond}*NModel.SGTSX{cCond};
    NModel.V21{cCond} = NModel.SGTSX{cCond}'*NModel.X21{cCond}*NModel.SFTSX;
end;

for (cCond1 = 1:N_CONDITIONS)
    for (cCond2 = 1:N_CONDITIONS)
        NModel.V22mDI{cCond1}{cCond2} = NModel.SGTSX{cCond1}'*NModel.X22mDI{cCond1}{cCond2}*NModel.SGTSX{cCond2};
    end;
end;

ScTerm = NModel.V11;

for (cCond = 1:N_CONDITIONS)
    ScTerm = ScTerm + NModel.V12{cCond};
    ScTerm = ScTerm + NModel.V21{cCond};
end;

for (cCond1 = 1:N_CONDITIONS)
    for (cCond2 = 1:N_CONDITIONS)
        ScTerm = ScTerm + NModel.V22mDI{cCond1}{cCond2};
    end;
end;

MhDist = NModel.SXTSX - ScTerm - NModel.SGXGI;

NModel.NTerm = NModel.NTerm - N_OBS_DIM*N_DATA*0.5*log(2*pi);
NModel.Loglike = NModel.NTerm - 0.5*NModel.LogFrDet - 0.5*(NModel.LogScDet1 + NModel.LogScDet2) - 0.5*MhDist;

