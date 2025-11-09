function [K_trim, row_idx, col_idx] = trim_nonzero(K)
% TRIM_NONZERO  Remove all-zero rows and columns from a matrix.
%
% Usage:
%   [K_trim, row_idx, col_idx] = trim_nonzero(K)
%
% Inputs:
%   K         - Input matrix (e.g., sparse or dense controller gain)
%
% Outputs:
%   K_trim    - Matrix with all-zero rows and columns removed
%   row_idx   - Indices of retained (nonzero) rows
%   col_idx   - Indices of retained (nonzero) columns
%
% Example:
%   K = [0 0 1; 0 0 0; 2 0 0];
%   [K_trim, r, c] = trim_nonzero(K);
%   % K_trim = [0 0 1; 2 0 0]
%   % r = [1 3], c = [1 3]
%
% Author: Pengyang Wu
% ================================================================

    % 找出非零行/列索引
    row_idx = find(any(K, 2));  % any nonzero across columns
    col_idx = find(any(K, 1));  % any nonzero across rows

    % 提取对应子矩阵
    K_trim = K(row_idx, col_idx);
end
