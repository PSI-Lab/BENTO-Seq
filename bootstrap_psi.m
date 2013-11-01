% estimate PSI by bootstrap
file_name = 'triplet_junc_rd.tsv'; tic
[ c1a,ac2,c1c2,id,n_c1a,n_ac2,n_c1c2,cover ] = load_TRI_RD_file(file_name);
fprintf('Loaded triplet read distribution file %s...\n', file_name);
n_inc = (n_c1a + n_ac2)/2; n_exc = n_c1c2;
psi_std = (n_inc + 1)./(n_inc + 1 + n_exc +1);
N = numel(id); K = 10000; %number of bootstrap samples
psi_bootstrap = zeros(N,1); psi_bootstrap_std = psi_bootstrap;
psi_samples = zeros(N,K);
for i = 1:N
    psit = gen_BootStrapSamples(c1a(i,:),ac2(i,:),c1c2(i,:),K);
    psi_samples(i,:) = psit;
    psi_bootstrap(i) = mean(psit);
    psi_bootstrap_std(i) = std(psit);
    if (mod(i,1000) == 0)
        fprintf('%d out of %d triplets processed...\n', i, N);
    end
end
toc