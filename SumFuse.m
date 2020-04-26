function [res, scores] = ProdFuse(frs, fri, ibas, idis, stinds, edinds)

N_SUBS = length(stinds); 
ND2 = size(frs{ibas}{idis}, 1);

scores = ones(ND2, N_SUBS)*-Inf; 
vdi = zeros(1, N_SUBS);

for csub = 1:N_SUBS
    cfrs = frs{ibas}{idis}(:, stinds(csub):edinds(csub));
    cfri = fri{ibas}(stinds(csub):edinds(csub));

    mxr = -ones(ND2, 1)*Inf;

    for i = 1:ND2
        mxr(i) = max(cfrs(i, :));
    end;

    nc = size(cfrs, 2);
    rpm = repmat(mxr, 1, nc);

    normfrs = cfrs - rpm;
    sumfrs = sum(exp(normfrs), 2);
    sumfrs = log(sumfrs) + mxr;

    scores(:, csub) = sumfrs;
    vdi(csub) = cfri(1);
end;

[mm pp] = max(scores);
tst = (pp == vdi);

res = sum(tst)/N_SUBS;

