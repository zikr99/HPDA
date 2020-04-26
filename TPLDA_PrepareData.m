function PData = TPLDA_PrepareData(Data, Conditions, Model, ID)

N_DATA = size(Data, 2);

for i = 1:N_DATA
    PData.XM(:, i) = Data(:, i) - Model.meanVec{Conditions(i)};

    SX = PData.XM(:, i).*Model.invSigmaEst{Conditions(i)}; 

    PData.FTSX(:, i) = Model.F{Conditions(i)}'*SX;
    PData.GTSX(:, i) = Model.G{ID}{Conditions(i)}'*SX;
    PData.XTSX(i) = PData.XM(:, i)'*SX;

    PData.GXGI(i) = PData.GTSX(:, i)'*Model.InvGTSGPlusI{ID}{Conditions(i)}*PData.GTSX(:, i);
end;
