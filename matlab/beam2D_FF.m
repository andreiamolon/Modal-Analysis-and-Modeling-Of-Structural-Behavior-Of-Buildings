% portique Locie
%----------------------------------------------------------------
% PURPOSE 
%    Set up the fe-model and perform eigenvalue analysis
%    for a simple frame structure.
%----------------------------------------------------------------

% portique Locie - liaisons semi-rigides - on crée des poutres dont les extrémités sont au niveau des accéléros

%----------------------------------------------------------------
% echo on;
 clear all;

% ------ Generate the model ------------------------------------------

% ------ material data ------------------------------------------
E=2e10;            rho=2500;
A=40e-3;           I=133e-6;            
nelem = 8;         ndof = 27;                 
% Boundary conditions
bc=[1 2 3 25 26 27]';

% ------ topology -----------------------------------------------
Edof=[1   1  2  3  4  5  6    
      2   4  5  6  7  8  9
      3   7  8  9  10  11  12    
      4   10 11 12 13  14  15
      5   13  14  15  16 17 18    
      6   16 17 18 19 20 21
      7   19 20 21 22 23 24   
      8   22 23 24 25 26 27
	  ];
% ------ list of coordinates  -----------------------------------
Coord=[0  0  ; 0.5 0 ; 1  0  ; 1.5 0 ; 2 0  ; 
    2.5 0  ; 3 0  ; 3.5 0 ; 4 0 
			];
% ------ list of degrees-of-fredom  -----------------------------
Dof=[1  2  3;  4  5  6; 7 8 9; 10 11 12; 13 14 15;16 17 18; 19 20 21; 22 23 24; 25 26 27];
% coordxtr = Extract element coordinates from a global coordinate matrix
 [Ex,Ey]=coordxtr(Edof,Coord,Dof,2);

% fixity factors
% : Create vectors and do matrix subscripting
% j: i: k is the same as [j,j + i, j + 2i, ... , k]
% de 0.6 ate 0.8 com um intervalo de 0.02
lp3=0.4:0.0001:1; 
%lp2=0:0.05:1;

err=1.e+6;
f1obj=35.207;
f2obj=99.858;
f3obj=180.219;

 i1=0;
 for p3=lp3
    i1=i1+1; 
    %i2=0;
    %for p2=lp2
        %i2=i2+1;
        
            % element properties 
            ep=[	E A I rho*A 1.0 1.0;
                    E A I rho*A 1.0 1.0;
                    E A I rho*A 1.0 1.0;
                    E A I rho*A 1.0 p3; 
                    E A I rho*A p3 1.0;
                    E A I rho*A 1.0 1.0;
                    E A I rho*A 1.0 1.0;
                    E A I rho*A 1.0 1.0;
                    ];
    % ------ generate element matrices, assemble in global matrices - 
              K=zeros(ndof);     M=zeros(ndof);
              
            for i=1:nelem
              [k,m,c]=beam2d_sr(Ex(i,:),Ey(i,:),ep(i,:));
              K=assem(Edof(i,:),K,k);  M=assem(Edof(i,:),M,m);  
            end

    % ----- Eigenvalue analysis --------------------------------------
    % Eigen = Solve a generalized eigenvalue problem
            [La,Egv]=eigen(K,M,bc);
             
            % formula para calcular a frequencia = sqrt(Lambida)/2*PI
            Freq1=sqrt(La(1))/(2*pi);
            disp(Freq1);
            Freq2=sqrt(La(2))/(2*pi);
            disp(Freq2);
            Freq3=sqrt(La(3))/(2*pi);
            disp(Freq3);
            % (f3obj-Freq3)^2
            err1=sqrt((f1obj-Freq1)^2+(f2obj-Freq2)^2+(f3obj-Freq3)^2)/3;
            if (err1<err)
                f1opt = Freq1; f2opt=Freq2; f3opt=Freq3; 
                p3opt=p3;%p2opt=p2;
                err=err1;
            end
        
    %end
 end

 % ----- plot one to 3 eigenmode ---------------------------------------
Freqall = [Freq1;Freq2;Freq3];
 for i=1:3
	figure(i),    clf,     grid,     title(strcat('Eigenmode ',int2str(i))), 
	eldraw2(Ex,Ey,[2 3 1]); 
	Edb=extract(Edof,Egv(:,i));      eldisp2(Ex,Ey,Edb,[1 2 2]);
	FreqText=num2str(Freqall(i));       text(.5,1.75,FreqText);
end
% -------------------- end ---------------------------------------
echo off;
