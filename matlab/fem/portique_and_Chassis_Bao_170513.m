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
E=2e+10;                 rho=2500;          %Betob
Av=50e-3;           Iv=260e-6;             % BA - 0.20x0.25
Ah=40e-3;           Ih=133e-6;             % BA - 0.20x0.20

EA=2e+11;                 rhoA=7850;          %Acier
A1=197.78e-4;           I1=57680.5e-8;             % Iy HEB400
A2=179.879e-4;           I2=40700e-8;             % Iy UUPN400
% Boundary conditions

bc=[1 2 34 35 88 89 90 91 92 93]';   %encastrements pour les noeuds 30,31-chassis
% bc=[1 2 3 34 35 36 88 89 90 91 92 93]';   %encastrements pour les noeuds 1, 12, 30 et 31-chassis
% bc=[1 2 3 34 35 36 46 47 48 76 77 78 88 89 90 91 92 93]';   %encastrements pour les noeuds 1, 12 16, 26, 30 et 31-chassis
% bc=[88 89 91 92]';   %rotules pour les noeuds 30 et 31 - chassis

%% ------ topology -----------------------------------------------
% numéro de barres et numéros des DDL des 2 extrémités
Edof=[1   1  2  3  4  5  6
      2   4  5  6  7  8  9
      3   7  8  9 10 11 12
      4  4 5 6 13 14 15
	  5  13 14 15 16 17 18
	  6  16 17 18 19 20 21
	  7  19 20 21 22 23 24
	  8  22 23 24 25 26 27
	  9  25 26 27 28 29 30
	  10 28 29 30 31 32 33
	  11 31 32 33 37 38 39
	  12 34 35 36 37 38 39
	  13 37 38 39 40 41 42
	  14 40 41 42 43 44 45
      15 46	47	48	49	50	51
      16	49	50	51	52	53	54
     17	52	53	54	55	56	57
    18	52	53	54	10	11	12
   19	10	11	12	58	59	60
   20	58	59	60	64	65	66
   21	64	65	66	70	71	72
   22	70	71	72	43	44	45
   23	43	44	45	82	83	84
   24	76	77	78	79	80	81
   25	79	80	81	82	83	84
   26	82	83	84	85	86	87
   27	1	2	3	61	62	63
   28	61	62	63	67	68	69
   29	67	68	69	73	74	75
   30	73	74	75	34	35	36
   31	46	47	48	88	89	90
   32	88	89	90	1	2	3
   33	34	35	36	91	92	93
   34	91	92	93	76	77	78
	  ];

% for n=1:31
%     Edof(n,:)=[n 3*n-2 3*n-1 3*n 3*n+1 3*n+2 3*n+3];
% end

%% ------ list of degrees-of-fredom  -----------------------------
% Dof=[1  2  3; 4  5  6; 7  8  9; 10 11 12; 13 14 15; 16 17 18; 19 20 21;
%      22 23 24; 25 26 27; 28 29 30; 31 32 33; 34 35 36; 
% 	 37 38 39; 40 41 42; 43 44 45];
for h=1:31
    Dof(h,:)=[3*(h-1)+1 3*(h-1)+2 3*(h-1)+3];
end
clear h;
% clear n h;

%% ------ list of coordinates  -----------------------------------
Coord=[0  0; 0  1.0; 0  1.55; 0  2; 0.465 1;
		0.645 1; 0.805 1; 1.26 1; 1.715 1; 1.875 1;
		2.055 1; 2.52 0; 2.52 1; 2.52 1.55; 2.52 2;
        -0.54 0; -0.54 1; -0.54 2; -0.54 2.6;   %noeuds 16 à 18
        0.63 2; 0.63 0;                        %noeuds 20 à 21
        1.26 2; 1.26 0;                        %noeuds 22 à 23
        1.89 2; 1.89 0;                        %noeuds 24 à 25
        3.06 0; 3.06 1; 3.06 2; 3.06 2.6        %noeuds 26 à 29
        -0.04 0; 2.56 0                         %noeuds 30 à 31
			];


%% fixity factors
%lp1=0.94; 
lp1=0.8; 
%lp1=0.8;
% lp1=0.7:0.02:0.98;
% lp2=0.4:0.02:1;
% lp2=0.665;
lp2=0.8;
% lp3=0.4:0.02:1;
lp3=0.84;

%lp4=0.4:0.02:1;
%lp4=0.5;
lp4=0.0;
p4=lp4;



% err=1.e+6;
%err=2;
%f1obj=73.24;f2obj= 85.45;f3obj=220.0;
%f1obj=70.19;f2obj= 85.00;f3obj=217.9;
%f1obj=70.19;f2obj= 85.00;f3obj=194.6;
%f1obj=61;f2obj= 85.00;f3obj=174.0;
% f1obj=73.24; f2obj= 85.45; f3obj=190.3;

p1=lp1;
p2=lp2;
p3=lp3;

          
            ep=[	E Av Iv rho*Av p1 1.0;                    %barre 1
                    E Av Iv rho*Av 1.0 1.0;                   %barre 2
                    E Av Iv rho*Av 1.0 p4;             %barre 3
                    E Ah Ih rho*Ah p2 1.0;              %barre 4
                    E Ah Ih rho*Ah 1.0 1.0;
                    E Ah Ih rho*Ah 1.0 1.0;
                    E Ah Ih rho*Ah 1.0 p3;
                    E Ah Ih rho*Ah 1.0 1.0;
                    E Ah Ih rho*Ah 1.0  1.0;
                    E Ah Ih rho*Ah 1.0  1.0;
                    E Ah Ih rho*Ah 1.0 p2;
                    E Av Iv rho*Av p1 1.0;
                    E Av Iv rho*Av 1.0 1.0;
                    E Av Iv rho*Av 1.0 p4;              %barre 14
                    EA A1 I1 rhoA*A1 1.0 1.0;    %barre 15
                    EA A1 I1 rhoA*A1 1.0 1.0;    %barre 16
                    EA A1 I1 rhoA*A1 1.0 1.0;    %barre 17
                    EA A1 I1 rhoA*A1 1.0 1.0;    %barre 18
                    EA A1 I1 rhoA*A1 1.0 1.0;    %barre 19
                    EA A1 I1 rhoA*A1 1.0 1.0;    %barre 20
                    EA A1 I1 rhoA*A1 1.0 1.0;    %barre 21
                    EA A1 I1 rhoA*A1 1.0 1.0;    %barre 22
                    EA A1 I1 rhoA*A1 1.0 1.0;    %barre 23
                    EA A2 I2 rhoA*A2 1.0 1.0;    %barre 24
                    EA A2 I2 rhoA*A2 1.0 1.0;    %barre 25
                    EA A2 I2 rhoA*A2 1.0 1.0;    %barre 26
                    EA A1 I1 rhoA*A1 1.0 1.0;    %barre 27
                    EA A1 I1 rhoA*A1 1.0 1.0;    %barre 28
                    EA A1 I1 rhoA*A1 1.0 1.0;    %barre 29
                    EA A1 I1 rhoA*A1 1.0 1.0;    %barre 30
                    EA A1 I1 rhoA*A1 1.0 1.0;    %barre 31
                    EA A1 I1 rhoA*A1 1.0 1.0;    %barre 32
                    EA A1 I1 rhoA*A1 1.0 1.0;    %barre 33
                    EA A1 I1 rhoA*A1 1.0 1.0;    %barre 34                    
                    ];
   
                %%
 [Ex,Ey]=coordxtr(Edof,Coord,Dof,2);
 
  
%%  ----- Draw a plot of the element mesh --------------------------


clf;     eldraw2(Ex,Ey,[1 2 2],Edof);
%grid;    
title('2-D Frame Structure') 
 
 %% ------ generate element matrices, assemble in global matrices - 
            K=zeros(93);     M=zeros(93);
            for i=1:31                                  %nombre des barres
              [k,m,c]=beam2d_sr(Ex(i,:),Ey(i,:),ep(i,:));
              K=assem(Edof(i,:),K,k);  M=assem(Edof(i,:),M,m);  
            end

    % ----- Eigenvalue analysis --------------------------------------

%             [La,Egv]=eigen(K,M,bc);
%             Freq1=sqrt(La(1))/(2*pi);
%             Freq2=sqrt(La(2))/(2*pi);
%             Freq3=sqrt(La(3))/(2*pi);
%              Freq4=sqrt(La(4))/(2*pi);
%               Freq5=sqrt(La(5))/(2*pi);
%                Freq6=sqrt(La(6))/(2*pi);
%                 Freq7=sqrt(La(7))/(2*pi);

[La,Egv]=eigen(K,M,bc);
Freq=sqrt(La)/(2*pi);

%% ----- plot eight eigenmodes ------------------------------------

figure(2), clf, 
% axis('equal'), 
axis('normal'), 
hold on, 
axis off
sfac=2;     %scale factor for displacements
title('The first eight eigenmodes (Hz)')
for i=1:4;
 Edb=extract(Edof,Egv(:,i));    %Edb: element displacement matrix
 Ext=Ex+(i-1)*3;                eldraw2(Ext,Ey,[1 2 1]); % draw the undeformed
 eldisp2(Ext,Ey,Edb,[2 3 2],sfac);  % draw the deformed 2D mesh
 FreqText=num2str(Freq(i));     text(3*(i-1)+.5,1.5,FreqText);
end;
Eyt=Ey-4; 
for i=5:8;
 Edb=extract(Edof,Egv(:,i));
 Ext=Ex+(i-5)*3;                eldraw2(Ext,Eyt,[1 2 1]); 
 eldisp2(Ext,Eyt,Edb,[2 3 2],sfac);
 FreqText=num2str(Freq(i));     text(3*(i-5)+.5,-2.5,FreqText);
end;


%%
%  i1=0;
%  for p1=lp1
%     i1=i1+1; 
%     i2=0;
%     for p2=lp2
%         i2=i2+1;
%         i3=0;
%         for p3=lp3
%             i3=i3+1;
%             
% %             p4=p1;
%             
%             ep=[	E Av Iv rho*Av p1 1.0;                    %barre 1
%                     E Av Iv rho*Av 1.0 1.0;                   %barre 2
%                     E Av Iv rho*Av 1.0 p4;             %barre 3
%                     E Ah Ih rho*Ah p2 1.0;              %barre 4
%                     E Ah Ih rho*Ah 1.0 1.0;
%                     E Ah Ih rho*Ah 1.0 1.0;
%                     E Ah Ih rho*Ah 1.0 p3;
%                     E Ah Ih rho*Ah 1.0 1.0;
%                     E Ah Ih rho*Ah 1.0  1.0;
%                     E Ah Ih rho*Ah 1.0  1.0;
%                     E Ah Ih rho*Ah 1.0 p2;
%                     E Av Iv rho*Av p1 1.0;
%                     E Av Iv rho*Av 1.0 1.0;
%                     E Av Iv rho*Av 1.0 p4;              %barre 14
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 15
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 16
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 17
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 18
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 19
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 20
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 21
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 22
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 23
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 24
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 25
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 26
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 27
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 28
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 29
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 30
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 31
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 32
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 33
%                     EA A1 I1 rhoA*A1 1.0 1.0;    %barre 34                    
%                     ];
%     % ------ generate element matrices, assemble in global matrices - 
%             K=zeros(93);     M=zeros(93);
%             for i=1:31                                  %nombre des barres
%               [k,m,c]=beam2d_sr(Ex(i,:),Ey(i,:),ep(i,:));
%               K=assem(Edof(i,:),K,k);  M=assem(Edof(i,:),M,m);  
%             end
% 
%     % ----- Eigenvalue analysis --------------------------------------
% 
%             [La,Egv]=eigen(K,M,bc);
%             Freq1=sqrt(La(1))/(2*pi);
%             Freq2=sqrt(La(2))/(2*pi);
%             Freq3=sqrt(La(3))/(2*pi);
%              Freq4=sqrt(La(4))/(2*pi);
% %               Freq5=sqrt(La(5))/(2*pi);
% %                Freq6=sqrt(La(6))/(2*pi);
% %                 Freq7=sqrt(La(7))/(2*pi);
%             err1=sqrt(9*(f1obj-Freq1)^2+4*(f2obj-Freq2)^2+(f3obj-Freq3)^2)/14;
%             if (err1<err)
%                 f1opt = Freq1; f2opt=Freq2 ;   f3opt=Freq3;
%                 p1opt=p1;p2opt=p2;p3opt=p3;
%                 err=err1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
%             end
%         end
%     end
%  end
 
%  %% ----- plot one to 3 eigenmode ---------------------------------------
% for i=1:3;
% 	figure(i),    clf,     grid,     title(strcat('Eigenmode ',int2str(i))), 
% 	eldraw2(Ex,Ey,[2 3 1]); 
% 	Edb=extract(Edof,Egv(:,i));      eldisp2(Ex,Ey,Edb,[1 2 2]);
% 	FreqText=num2str('f',i,'opt'));       text(.5,1.75,FreqText);
% end
% -------------------- end ---------------------------------------
echo off;
