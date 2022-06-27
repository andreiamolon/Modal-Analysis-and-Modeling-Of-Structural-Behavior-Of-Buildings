% portique Locie
%----------------------------------------------------------------
% PURPOSE 
%    Set up the fe-model and perform eigenvalue analysis
%    for a simple frame structure.
%----------------------------------------------------------------

% portique Locie - liaisons semi-rigides 
% modèle simplifié 6"poutres" : 2 par poteau, 2 dans la poutre
%----------------------------------------------------------------
 echo on

% ------ Generate the model ------------------------------------------

% ------ material data ------------------------------------------
E=2e10;             
rho=2500;           A= 40e-3;             % BA - 0.20x0.25
m = rho*A;          I=133e-6;             % BA - 0.20x0.20
nelem = 4;          ndof = 15;
%rigidité des liaisons
p1=1; p2=1; %p1z=1.0; p2z=1.0;
% Ep = element properties 
%    E A I m  p1 p2  
ep=[	E A I m p1 p2;
        E A I m p1 p2; 
        E A I m p1  p2;
        E A I m p1 p2;
	];
% ------ topology -----------------------------------------------
Edof=[1   1  2  3  4  5  6    
      2   4  5  6  7  8  9
      3   7  8  9  10  11  12    
      4   10 11 12 13  14  15  
     ];
% ------ list of coordinates  -----------------------------------
Coord=[0  0 ; 1  0 ; 2 0 ; 3 0 ; 4 0
	   ];
% ------ list of degrees-of-fredom  -----------------------------
Dof=[1  2  3;  4  5  6; 7 8 9; 10 11 12; 13 14 15
      	];
% ------ generate element matrices, assemble in global matrices - 
K=zeros(ndof);     M=zeros(ndof);
[Ex,Ey] = coordxtr(Edof,Coord,Dof,2); 
% eo = [xz yz zz] orientation of local z axis

for i=1:nelem
  %eo(i,:) = [0 0 1];
  [k,m,c]=beam2d_sr(Ex(i,:),Ey(i,:),ep(i,:));
  K=assem(Edof(i,:),K,k);  M=assem(Edof(i,:),M,m);  

end
 
% ----- Draw a plot of the element mesh --------------------------

% Draw undeformed finite element mesh
% search for a function to draw mesh 3D
clf;     eldraw2(Ex,Ey,[2 3 1],Edof);
grid;    title('2-D Frame Structure') 

% ----- Eigenvalue analysis --------------------------------------

b=[1 2 3 13 14 15]';
[La,Egv]=eigen(K,M,b);
Freq=sqrt(La)/(2*pi);
Freqreal=real(Freq);
n_modes= length (Freq);

% ----- plot one to 3 eigenmode ---------------------------------------
%for i=1:3
%   	figure(i),    clf,     grid,     title(strcat('Eigenmode ',int2str(i))), 
%	eldraw2(Ex,Ey,[2 3 1]); 
%	Edb=extract(Edof,Egv(:,i));      eldisp2(Ex,Ey,Edb,[2 3 1]);
%	FreqText=num2str(Freq(i));       text(.5,1.75,FreqText);
%end
% ----- plot eight eigenmodes ------------------------------------

%figure(2), clf, axis('equal'), hold on, axis off
%sfac=0.5;
%title('The first eight eigenmodes (Hz)' )
%for i=1:4;
%  Edb=extract(Edof,Egv(:,i));
%  Ext=Ex+(i-1)*3;                eldraw2(Ext,Ey,[2 3 1]); 
%  eldisp2(Ext,Ey,Edb,[1 2 2],sfac);
%  FreqText=num2str(Freq(i));     text(3*(i-1)+.5,1.5,FreqText);
%end;
%Eyt=Ey-4; 
%for i=5:8;
%  Edb=extract(Edof,Egv(:,i));
%  Ext=Ex+(i-5)*3;                eldraw2(Ext,Eyt,[2 3 1]); 
%  eldisp2(Ext,Eyt,Edb,[1 2 2],sfac);
%  FreqText=num2str(Freq(i));     text(3*(i-5)+.5,-2.5,FreqText);
%end

% -------------------- end ---------------------------------------
echo off
