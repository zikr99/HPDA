function result = OnWorkingModel(WModel, PData, condition, ID)

N_CONDITIONS = length(WModel.F);
N_HID_DIM = size(WModel.F{1}, 2);
N_HID_DIM_NOISE = size(WModel.G{1}{1}, 2);
N_OBS_DIM = size(PData.XM, 1);
N_DATA = size(PData.XM, 2);

WModel.X00 = WModel.X00 + WModel.AiMinBDICi{ID}{condition};

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

WModel.LogScDet1 = logdet(WModel.X00);

WModel.LogFrDet = WModel.LogFrDet + sum(log(WModel.Sigma{condition}));
WModel.LogScDet2 = WModel.LogScDet2 + WModel.LogDetGTSGPlusI{ID}{condition};

WModel.SFTSX = repmat(WModel.SFTSX, 1, N_DATA) + PData.FTSX;

for (cCond = 1:N_CONDITIONS)
    WModel.SGTSX{cCond} = repmat(WModel.SGTSX{cCond}, 1, N_DATA);
end;

WModel.SGTSX{condition} = WModel.SGTSX{condition} + PData.GTSX;

WModel.SXTSX = repmat(WModel.SXTSX, 1, N_DATA) + PData.XTSX;
WModel.SGXGI = repmat(WModel.SGXGI, 1, N_DATA) + PData.GXGI;  

WModel.V11 = sum((WModel.SFTSX'*WModel.X11).*WModel.SFTSX', 2);

for (cCond = 1:N_CONDITIONS)
    WModel.V12{cCond} = sum((WModel.SFTSX'*WModel.X12{cCond}).*WModel.SGTSX{cCond}', 2);
    WModel.V21{cCond} = sum((WModel.SGTSX{cCond}'*WModel.X21{cCond}).*WModel.SFTSX', 2);
end;

for (cCond1 = 1:N_CONDITIONS)
    for (cCond2 = 1:N_CONDITIONS)
        WModel.V22mDI{cCond1}{cCond2} = sum((WModel.SGTSX{cCond1}'*WModel.X22mDI{cCond1}{cCond2}).*WModel.SGTSX{cCond2}', 2);
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

MhDist = WModel.SXTSX - ScTerm' - WModel.SGXGI;

WModel.NTerm = WModel.NTerm - N_OBS_DIM*1*0.5*log(2*pi);
WModel.Loglike = WModel.NTerm - 0.5*WModel.LogFrDet - 0.5*(WModel.LogScDet1 + WModel.LogScDet2) - 0.5*MhDist;

result = WModel.Loglike;

