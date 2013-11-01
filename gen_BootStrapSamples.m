function [ psit ] = gen_BootStrapSamples(c1a,ac2,c1c2,K)
n=size(c1a,2);
ninc=(sum(c1a(ceil(rand(n,K)*n)))+sum(ac2(ceil(rand(n,K)*n))));
nexc=(sum(c1c2(ceil(rand(n,K)*n))));
incl=randg(ninc/2+1);
excl=randg(nexc+1);
all_l=incl+excl;
psit=incl./all_l;
end
