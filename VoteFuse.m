function [res, scores] = VoteFuse(frs, fri, ibas, idis, stinds, edinds)

N_SUBS = length(stinds); 
ND2 = size(frs{ibas}{idis}, 1);

scores = zeros(ND2, N_SUBS); 
vdi = zeros(1, N_SUBS);

for csub = 1:N_SUBS
    cfrs = frs{ibas}{idis}(:, stinds(csub):edinds(csub));
    cfri = fri{ibas}(stinds(csub):edinds(csub));

    nc = edinds(csub) - stinds(csub) + 1;
    tpvote = zeros(ND2, nc); 

    [m p] = max(cfrs);
    
    for i = 1:nc
        tpvote(p(i), i) = 1;
    end; 

    scores(:, csub) = sum(tpvote, 2);;
    vdi(csub) = cfri(1);
end;

[mm pp] = max(scores);
tst = (pp == vdi);

res = sum(tst)/N_SUBS;

