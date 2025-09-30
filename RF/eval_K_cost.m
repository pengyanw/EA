function [cost, isstable,specrad] = eval_K_cost(K_sparse, A, B_, Q_bm, R_bm_, costBM, alpha_cost, max_links)
% alpha*lqrcost+(1-alpha)*nnz(K_sparse)/max_links
lqr_cost = get_lqr_cost(A,B_,Q_bm,R_bm_,K_sparse)/costBM;
sparse_cost = nnz(K_sparse)/max_links;
cost = alpha_cost*lqr_cost + (1-alpha_cost)*sparse_cost;
isstable = lqr_cost>=1e2;
specrad = max(abs(eig(A+B_*K_sparse)));