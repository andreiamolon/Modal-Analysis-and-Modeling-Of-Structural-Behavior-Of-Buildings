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
ma = rho*A;         Iz=133e-6;            nelem = 8;         ndof = 54;           
% Boundary conditions
bc=[1 2 3 4 5 6 10 16 22 28 34 40 46 49 50 51 52 53 54]';

% ------ topology -----------------------------------------------
Edof=[1   1  2  3  4  5  6  7  8  9 10 11 12  
      2   7  8  9 10 11 12 13 14 15 16 17 18
      3   13 14 15 16 17 18 19 20 21 22 23 24
      4   19 20 21 22 23 24 25 26 27 28 29 30
      5   25 26 27 28 29 30 31 32 33 34 35 36
      6   31 32 33 34 35 36 37 38 39 40 41 42
      7   37 38 39 40 41 42 43 44 45 46 47 48
      8   43 44 45 46 47 48 49 50 51 52 53 54
     ];
% ------ list of coordinates  -----------------------------------
Coord=[0  0  0 ; 0.25 0 0; 0.5  0  0 ; 0.75 0 0; 1 0 0 ; 
    1.25 0 0 ; 1.5 0 0 ; 1.75 0 0; 2 0 0; 
	   ];
% ------ list of degrees-of-fredom  -----------------------------
Dof=[1  2  3  4  5  6; 7 8 9 10 11 12; 13 14 15 16 17 18 ; 19 20 21 22 23 24; 25 26 27 28 29 30;
      	31 32 33 34 35 36; 37 38 39 40 41 42; 43 44 45 46 47 48; 49 50 51 52 53 54];
[Ex,Ey,Ez] = coordxtr(Edof,Coord,Dof,2);

% fixity factors
% : Create vectors and do matrix subscripting
% j: i: k is the same as [j,j + i, j + 2i, ... , k]
% de 0.6 ate 0.8 com um intervalo de 0.02
%lp1y=0:0.05:1;
%lp1z=0:0.05:1;
%lp2y=0:0.05:1;
%lp2z=0:0.05:1;
lp3y=0.5:0.001:1;
lp3z=0.5:0.001:1;

err=1.e+6;
%freq theorique 

f1obj = 139.35;
f2obj = 139.35;
f3obj = 359.941;
f4obj = 359.941;
f5obj = 655.813;
f6obj= 655.813;


i1=0;
 for p3y=lp3y
    i1=i1+1; 
    i2=0;
    lp3z=0.5:0.001:p3y;
    for p3z=lp3z
        i2=i2+1;
        %i3=0;
       % for p1z=lp1z
            %i3=i3+1;
            %i4=0;
           % for p2z=lp2z
               % i4=i4+1;
% 
                % element properties 
                ep=[ E G A Iy Iz ma Kv 1.0 1.0 1.0 1.0;
                     E G A Iy Iz ma Kv 1.0 1.0 1.0 1.0;
                     E G A Iy Iz ma Kv 1.0 1.0 1.0 1.0;
		             E G A Iy Iz ma Kv 1.0 p3y 1.0 p3z;
                     E G A Iy Iz ma Kv p3y 1.0 p3z 1.0;
                     E G A Iy Iz ma Kv 1.0 1.0 1.0 1.0;
                     E G A Iy Iz ma Kv 1.0 1.0 1.0 1.0;
                     E G A Iy Iz ma Kv 1.0 1.0 1.0 1.0 ;
		                                     
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
    
        % ----- Eigenvalue analysis --------------------------------------
        % Eigen = Solve a generalized eigenvalue problem
        % Calculamos a frequencia e comparamos com a frequencia obtida em
        % laboratorio para calcular os fixity factors.
                [La,Egv]=eigen(K,M,bc);
                % formula para calcular a frequencia = sqrt(Lambida)/2*PI
                Freq1=sqrt(La(1))/(2*pi);                
                Freq1real = real(Freq1);
                Freq2=sqrt(La(2))/(2*pi);               
                Freq2real = real(Freq2);
                Freq3=sqrt(La(3))/(2*pi);                
                Freq3real = real(Freq3);
                Freq4=sqrt(La(4))/(2*pi);                
                Freq4real = real(Freq4);
                Freq5=sqrt(La(5))/(2*pi);                
                Freq5real = real(Freq5);
                Freq6=sqrt(La(6))/(2*pi);                
                Freq6real = real(Freq6);

                err1=(sqrt(((f1obj-Freq1)^2) / f1obj^2 +((f2obj-Freq2)^2) / f2obj^2 ))/2;
                if (err1<err)
                    f1opt=Freq1; 
                    f2opt=Freq2;
                    f3opt=Freq3;  
                    f4opt=Freq4; 
                    f5opt=Freq5; 
                    f6opt=Freq6;    
                    %p1yopt=p1y;p2yopt=p2y;p1zopt=p1z;p2zopt=p2z; 
                    p3yopt=p3y;p3zopt=p3z;
                    err=err1;
                end
            %end
        %end
    end
 end
 
% ----- Draw a plot of the element mesh --------------------------

% Draw undeformed finite element mesh
% search for a function to draw mesh 3D
clf;     eldraw3(Ex,Ey,Ez,[2 3 1],Edof);
grid;    title('3-D Frame Structure') 

% ----- plot one to 3 eigenmode ---------------------------------------
Freq = [Freq1;Freq2;Freq3;Freq4;Freq5;Freq6];
indexes = [1;2;3;4;5;6];
for i=1:3    
	figure(i),    clf,     grid,     title(strcat('Eigenmode ',int2str(i))), 
	eldraw3(Ex,Ey,Ez,[2 3 1]); 
	Edb=extract(Edof,Egv(:,indexes(i)));      eldisp3(Ex,Ey,Ez,Edb,[2 3 1]);
	FreqText=num2str(Freq(i));       text(.5,1.75,FreqText);
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
