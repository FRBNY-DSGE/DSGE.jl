%%% Code used to generate mat files read in in test_Gensys.jl

clear()

% G0, G1, etc. matrices from evaluating code up to gensys call 
load('gensys_args.mat')
[G1_gensys, C_gensys, impact, fmat, fwt, ywt, gev, eu] = gensys(G0, G1, C, PSI, PIE, div);

[AA, BB, Q, Z] = qz(G0, G1);
alpha = diag(AA);
beta = diag(BB);

[AA_qzdiv, BB_qzdiv, Q_qzdiv, Z_qzdiv] = qzdiv(div, AA, BB, Q, Z);

E = ordeig(AA, BB);
select = abs(E) < div;
[AA_ordqz, BB_ordqz, Q_ordqz, Z_ordqz] = ordqz(AA, BB, Q, Z, select);
%[AAS, BBS, QS, ZS] = ordqz(AA, BB, Q, Z, 'udi');

save('gensys_variables.mat')