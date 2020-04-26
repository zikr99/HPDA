function WModel = InitWorkingModel(PModel, PInitData, PInitConditions, ID)

N_CONDITIONS = length(PModel.F);
N_HID_DIM = size(PModel.F{1}, 2);
N_HID_DIM_NOISE = size(PModel.G{1}{1}, 2);
N_OBS_DIM = size(PInitData.XM, 1);
N_DATA = size(PInitData.XM, 2);

WModel = PModel;

WModel.X00 = zeros(N_HID_DIM, N_HID_DIM);

for (i = 1:N_DATA)
    WModel.X00 = WModel.X00 + WModel.AiMinBDICi{ID}{PInitConditions(i)};
end;

WModel.X11 = inv(WModel.X00 + eye(N_HID_DIM, N_HID_DIM));

for (cCond = 1:N_CONDITIONS)
    WModel.X12{cCond} = -WModel.X11*WModel.BDI{ID}{cCond};
    WModel.X21{cCond} = -WModel.DIC{ID}{cCond}*WModel.X11;
end;

for (cCond1 = 1:N_CONDITIONS)
    for (cCond2 = 1:N_CONDITIONS)
        WModel.X22mDI{cCond1}{cCond2} = -WModel.X21{cCond1}*WModel.BDI{ID}{cCond2};
    end;
end;

WModel.LogFrDet = 0;
WModel.LogScDet1 = logdet(WModel.X00);
WModel.LogScDet2 = 0;

for (i = 1:N_DATA)
    WModel.LogFrDet = WModel.LogFrDet + sum(log(WModel.Sigma{PInitConditions(i)}));
    WModel.LogScDet2 = WModel.LogScDet2 + WModel.LogDetGTSGPlusI{ID}{PInitConditions(i)};
end;

WModel.SFTSX = zeros(N_HID_DIM, 1);

for (i = 1:N_DATA)
    WModel.SFTSX = WModel.SFTSX + PInitData.FTSX(:, i);
end;

for (cCond = 1:N_CONDITIONS)
    WModel.SGTSX{cCond} = zeros(N_HID_DIM_NOISE, 1);

    iIndices = find(PInitConditions == cCond);

    for (i = iIndices')
        WModel.SGTSX{cCond} = WModel.SGTSX{cCond} + PInitData.GTSX(:, i);  
    end;
end;

WModel.SXTSX = 0;
WModel.SGXGI = 0;

for (i = 1:N_DATA)
    WModel.SXTSX = WModel.SXTSX + PInitData.XTSX(i);
    WModel.SGXGI = WModel.SGXGI + PInitData.GXGI(i);
end;

WModel.V11 = WModel.SFTSX'*WModel.X11*WModel.SFTSX;

for (cCond = 1:N_CONDITIONS)
    WModel.V12{cCond} = WModel.SFTSX'*WModel.X12{cCond}*WModel.SGTSX{cCond};
    WModel.V21{cCond} = WModel.SGTSX{cCond}'*WModel.X21{cCond}*WModel.SFTSX;
end;

for (cCond1 = 1:N_CONDITIONS)
    for (cCond2 = 1:N_CONDITIONS)
        WModel.V22mDI{cCond1}{cCond2} = WModel.SGTSX{cCond1}'*WModel.X22mDI{cCond1}{cCond2}*WModel.SGTSX{cCond2};
    end;
end;

ScTerm = WModel.V11;

for (cCond = 1:N_CONDITIONS)
    ScTerm = ScTerm + WModel.V12{cCond};
    ScTerm = ScTerm + WModel.V21{cCond};
end;

for (cCond1 = 1:N_CONDITIONS)
    for (cCond2 = 1:N_CONDITIONS)
        ScTerm = ScTerm + WModel.V22mDI{cCond1}{cCond2};
    end;
end;

MhDist = WModel.SXTSX - ScTerm - WModel.SGXGI;

WModel.NTerm = -N_OBS_DIM*N_DATA*0.5*log(2*pi);
WModel.Loglike = WModel.NTerm - 0.5*WModel.LogFrDet - 0.5*(WModel.LogScDet1 + WModel.LogScDet2) - 0.5*MhDist;

