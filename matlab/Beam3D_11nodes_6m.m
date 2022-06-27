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
E=2e10;         G = E/(2*(1+0.2));      
rho=2500;       Kv=0; %Saint-Venants torsion               
nelem = 10;
ndf = 66;
% Geometry 
A= 40e-3;       % A - 0.20x0.20
m = rho*A;
% Element Properties
Iz=133e-6;      Iy=133e-6;
%E G A Iy Iz Kv m  p1y p2y p1z p2z 
%rigidité des liaisons
p3y = 0.859; p3z = 0.787;
ep=[E G A Iy Iz ma Kv 1.0 1.0 1.0 1.0;
                     E G A Iy Iz ma Kv 1.0 1.0 1.0 1.0;
                     E G A Iy Iz ma Kv 1.0 1.0 1.0 1.0;
                     E G A Iy Iz ma Kv 1.0 1.0 1.0 1.0;
		             E G A Iy Iz ma Kv 1.0 p3y 1.0 p3z;
                     E G A Iy Iz ma Kv p3y 1.0 p3z 1.0;
                     E G A Iy Iz ma Kv 1.0 1.0 1.0 1.0;
                     E G A Iy Iz ma Kv 1.0 1.0 1.0 1.0;
                     E G A Iy Iz ma Kv 1.0 1.0 1.0 1.0;
                     E G A Iy Iz ma Kv 1.0 1.0 1.0 1.0 ;
	];
% ------ topology -----------------------------------------------
Edof=[1   1  2  3  4  5  6  7  8  9 10 11 12  
      2   7  8  9 10 11 12 13 14 15 16 17 18
      3   13 14 15 16 17 18 19 20 21 22 23 24
      4   19 20 21 22 23 24 25 26 27 28 29 30
      5   25 26 27 28 29 30 31 32 33 34 35 36
      6   31 32 33 34 35 36 37 38 39 40 41 42
      7   37 38 39 40 41 42 43 44 45 46 47 48
      8   43 44 45 46 47 48 49 50 51 52 53 54
      9   49 50 51 52 53 54 55 56 57 58 59 60
      10  55 56 57 58 59 60 61 62 63 64 65 66 
     ];
% ------ list of coordinates  -----------------------------------
Coord=[0  0  0 ; 0.8 0 0; 1.6  0  0 ; 2.4 0 0; 3.2 0 0 ; 
    4 0 0 ; 4.8 0 0 ; 5.6 0 0; 6.4 0 0; 7.2 0 0; 8 0 0]; 

% ------ list of degrees-of-fredom  -----------------------------
Dof=[1  2  3  4  5  6; 7 8 9 10 11 12; 13 14 15 16 17 18 ; 19 20 21 22 23 24; 25 26 27 28 29 30;31 32 33 34 35 36;37 38 39 40 41 42;43 44 45 46 47 48; 49 50 51 52 53 54;55 56 57 58 59 60; 61 62 63 64 65 66
      	];
% ------ generate element matrices, assemble in global matrices - 
K=zeros(ndf);     M=zeros(ndf);
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

b=[1 2 3 4 5 6 10 16 22 28 34 40 46 52 58 61 62 63 64 65 66]';
[L,Egv]=eigen(K,M,b);
Freq=sqrt(L)/(2*pi);
Freq1 = Freq(1);
Freq2 = Freq(2);
Freq3 = Freq(3);
Freq4 = Freq(4);

f1obj =35.514;
f2obj =35.813;
f3obj =97.150;
f4obj =97.150;
err1=(sqrt(((f1obj-Freq1)^2) / f1obj^2 +((f2obj-Freq2)^2) / f2obj^2 +((f3obj-Freq3)^2) / f3obj^2 + ((f4obj-Freq4)^2) / f4obj^2 ))/4;
                
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
