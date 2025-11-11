function ICC_mat = calculate_icc_matrix_three_sessions(data1, data2, data3)
    % Calculate ICC(3,1) for each edge across three sessions
    % data1, data2, data3: [node x node x subject]
    [n1, n2, nsubj] = size(data1);
    ICC_mat = nan(n1, n2);
    for i = 1:n1
        for j = 1:n2
            x = squeeze([data1(i,j,:); data2(i,j,:); data3(i,j,:)])'; % [nsubj x 3]
            if all(~isnan(x(:)))
                ICC_mat(i,j) = ICC_shrout3(x);
            end
        end
    end
end