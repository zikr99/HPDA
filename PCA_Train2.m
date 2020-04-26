function [PrComp SSq] = PCA_Train2(data, numcomp)

N = size(data, 2);
DIM = size(data, 1);

if (DIM < N)
    XTX = (1/N)*(data*data');
else
    XTX = (1/N)*(data'*data);
end

[dummy LSq V] = svd(XTX);
PrComp = zeros(size(data, 1), numcomp);

for (c = 1:numcomp)
    if (DIM < N)
        PrComp(:, c) = V(:, c);
    else
        PrComp(:, c) = data*V(:, c);
    end;

    PrComp(:, c) = PrComp(:, c)*(1/norm(PrComp(:, c)));
end;

EV = diag(LSq);
EV = EV(1:numcomp);

SSq = (trace(XTX) - sum(EV))/(size(data, 1) - numcomp);

for (c = 1:numcomp)
    L(c) = sqrt(EV(c) - SSq);
end;

for (c = 1:numcomp)
    PrComp(:, c) = PrComp(:, c)*L(c);
end;
