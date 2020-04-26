function result = OnSimpleModel(SPModel, PData, condition, ID)

N_CONDITIONS = length(SPModel.F);
N_HID_DIM = size(SPModel.F{1}, 2);
N_HID_DIM_NOISE = size(SPModel.G{1}{1}, 2);
N_OBS_DIM = size(PData.XM, 1);
N_DATA = size(PData.XM, 2);

ATSX = [PData.FTSX; PData.GTSX];
ScTerm = sum((ATSX'*SPModel.InvIATSA{ID}{condition}).*ATSX', 2);

MhDist = PData.XTSX - ScTerm';

result = SPModel.FrLogTerm{ID}{condition} - 0.5*MhDist;

 