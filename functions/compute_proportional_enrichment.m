function [propTbl, statTbl] = compute_proportional_enrichment(topMask, netIndex)
% topMask: logical 90x90 for top variable edges (symmetric, lower used)
% netIndex: 90x1 with labels 1..7
    nNet = max(netIndex);
    maskL = tril(true(size(topMask,1)), -1);
    topL = false(size(topMask)); topL(maskL) = topMask(maskL);

    % Count possible edges per (i,j) network pair and observed top edges
    Obs = zeros(nNet,nNet);
    Pos = zeros(nNet,nNet);
    for a = 1:nNet
        for b = a:nNet
            idxA = find(netIndex==a);
            idxB = find(netIndex==b);
            if a==b
                posEdges = nchoosek(numel(idxA),2);
                obsEdges = nnz(topL(idxA, idxA));
            else
                posEdges = numel(idxA)*numel(idxB);
                obsEdges = nnz(topL(idxA, idxB));
            end
            Obs(a,b)=obsEdges; Obs(b,a)=obsEdges;
            Pos(a,b)=posEdges; Pos(b,a)=posEdges;
        end
    end
    Prop = Obs ./ max(Pos,1);
    propTbl = table(); % flatten to long form
    rows = {};
    for a = 1:nNet
        for b = a:nNet
            rows{end+1,1} = sprintf('%d-%d',a,b); %#ok<AGROW>
        end
    end
    propTbl.Pair = rows;
    propTbl.Observed = [];
    propTbl.Possible = [];
    propTbl.Proportion = [];
    k=1;
    for a=1:nNet
        for b=a:nNet
            propTbl.Observed(k,1) = Obs(a,b);
            propTbl.Possible(k,1) = Pos(a,b);
            propTbl.Proportion(k,1) = Prop(a,b);
            k=k+1;
        end
    end

    % Simple enrichment vs global rate using binomial test or hypergeometric
    totalObs = sum(Obs(tril(true(nNet))));
    totalPos = sum(Pos(tril(true(nNet))));
    global_rate = totalObs / totalPos;

    pvals = zeros(height(propTbl),1);
    k=1;
    for a=1:nNet
        for b=a:nNet
            x = Obs(a,b); m = Pos(a,b);
            % Binomial test: P[X >= x | m, p=global_rate]
            pvals(k) = 1 - binocdf(x-1, m, global_rate);
            k=k+1;
        end
    end
    [~, ~, q] = fdr_bh(pvals);
    statTbl = propTbl;
    statTbl.p_binom = pvals;
    statTbl.q_fdr  = q;
end