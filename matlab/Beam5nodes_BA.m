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
E=2e10;             G = E/2*(1+0.2);      rho=2500;           
A= 40e-3;           Iy=133e-6;            Kv=0;  %Saint-Venants torsion          
ma = rho*A;         Iz=133e-6;            nelem = 4;         ndof = 30; 
%rigidité des liaisons
p3y=1.0; p3z=1.0; 

 ep= [  E G A Iy Iz ma Kv 1.0 1.0 1.0 1.0;
        E G A Iy Iz ma Kv 1.0 p3y 1.0 p3z;
        E G A Iy Iz ma Kv p3y 1.0 p3z 1.0;
        E G A Iy Iz ma Kv 1.0 1.0 1.0 1.0;
		
	];
% ------ topology -----------------------------------------------
Edof=[1   1  2  3  4  5  6  7  8  9 10 11 12  
      2   7  8  9 10 11 12 13 14 15 16 17 18
      3   13 14 15 16 17 18 19 20 21 22 23 24
      4   19 20 21 22 23 24 25 26 27 28 29 30

     ];
% ------ list of coordinates  -----------------------------------
Coord=[0  0  0 ; 1 0 0 ; 2  0  0 ; 3 0 0; 4 0 0 ; 
   ]; 

% ------ list of degrees-of-fredom  -----------------------------
Dof=[1  2  3  4  5  6; 7 8 9 10 11 12; 13 14 15 16 17 18 ; 19 20 21 22 23 24; 25 26 27 28 29 30;
      	];
% ------ generate element matrices, assemble in global matrices - 
K=zeros(ndof);     M=zeros(ndof);
[Ex,Ey,Ez] = coordxtr(Edof,Coord,Dof,2); 
% eo = [xz yz zz] orientation of local z axis
    
    for i=1:nelem
      eo(i,:) = [0 0 1]; 
      [k,m,c]=beam3d_sr(Ex(i,:),Ey(i,:),Ez(i,:),eo(i,:),ep(i,:));
      K=assem(Edof(i,:),K,k);  M=assem(Edof(i,:),M,m);  
    
    end
 
% ----- Draw a plot of the element mesh --------------------------

% Draw undeformed finite element mesh
clf;     eldraw3(Ex,Ey,Ez,[2 3 1],Edof);
grid;    title('3-D Frame Structure') 

% ----- Eigenvalue analysis --------------------------------------

b=[1 2 3 4 10 16 22 28 25 26 27 ]';
[L,Egv]=eigen(K,M,b);
Freq=sqrt(L)/(2*pi);
freqreal = real(Freq);
n_modes= length (freqreal);

% ----- plot one to 3 eigenmode ---------------------------------------
for i=1:3
	figure(i),    clf,     grid,     title(strcat('Eigenmode ',int2str(i))), 
	eldraw3(Ex,Ey,Ez,[2,3,1]); 
	Edb=extract(Edof,Egv(:,i));
    eldisp3(Ex,Ey,Ez,Edb,[2 3 1]);
	FreqText=num2str(freqreal(i));       text(.5,1.75,FreqText);
end
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
