function [similarity_matrix] = calculate_similarity_matrix_seguin(data_session1, data_session2)
% Seguin 2022 node-wise correlation: for each node, correlate rows, then average across nodes.
[n_roi, ~, n_subj] = size(data_session1);
similarity_matrix = zeros(n_subj, n_subj);
for i = 1:n_subj
    A1 = squeeze(data_session1(:,:,i));
    for j = 1:n_subj
        A2 = squeeze(data_session2(:,:,j));
        row_corrs = zeros(n_roi,1);
        for k = 1:n_roi
            row_corrs(k) = corr(A1(k,:)', A2(k,:)');
        end
        similarity_matrix(i,j) = mean(row_corrs);
    end
end
end