function result = OnSimpleModel12(SPModel, PData1, PData2, condition, ID)

N_CONDITIONS = length(SPModel.F);
N_HID_DIM = size(SPModel.F{1}, 2);
N_HID_DIM_NOISE = size(SPModel.G{1}{1}, 2);
N_OBS_DIM = size(PData1.XM, 1);
N_DATA1 = size(PData1.XM, 2);
N_DATA2 = size(PData2.XM, 2);

result = -Inf*ones(N_DATA2, N_DATA1);

for c1 = 1:N_DATA1
    ATSX = [repmat(PData1.FTSX(:, c1), 1, N_DATA2) + PData2.FTSX; repmat(PData1.GTSX(:, c1), 1, N_DATA2); PData2.GTSX];
    ScTerm = sum((ATSX'*SPModel.InvIATSA{ID}{condition}).*ATSX', 2);

    MhDist = PData2.XTSX + PData1.XTSX(c1) - ScTerm';

    result(:, c1) = (SPModel.FrLogTerm{ID}{condition} - 0.5*MhDist)';
end;

