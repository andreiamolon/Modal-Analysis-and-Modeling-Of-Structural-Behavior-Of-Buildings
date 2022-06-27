
%----------------------------------------------------------------
% PURPOSE 
%    Set up the fe-model and perform eigenvalue analysis
%    for a simple frame structure.
%----------------------------------------------------------------

% REFERENCES
%     G"oran Sandberg 1994-03-08 
%     Karl-Gunnar Olsson 1995-09-29
%----------------------------------------------------------------

hauteur=0.5475;
dif=100;
erreur = 0.001;
f_castem=82.686;
mode=1;
while(dif > erreur && hauteur>0.5475)
    f = prd_hf(hauteur,mode);
    dif = abs(f-f_castem);
   hauteur=hauteur-0.00001;
   X=[f hauteur dif];
   disp(X) 
end

disp(hauteur)


