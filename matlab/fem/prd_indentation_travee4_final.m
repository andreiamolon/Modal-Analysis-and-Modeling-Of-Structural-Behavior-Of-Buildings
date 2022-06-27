
%----------------------------------------------------------------
% PURPOSE 
%    Set up the fe-model and perform eigenvalue analysis
%    for a simple frame structure.
%----------------------------------------------------------------

% REFERENCES
%     G"oran Sandberg 1994-03-08 
%     Karl-Gunnar Olsson 1995-09-29
%----------------------------------------------------------------
 echo off

% ------ Generate the model ------------------------------------------

% ------ material data ------------------------------------------
Eh=3e10;  Ev=2.7e10;    rho=2400;
Av=0.36e-2;           Iv=1.08e-6;             % Poteaux et poutres
Ah=0.0375;           Ih=1.953125e-4;            % Longrines

% ------ topology -----------------------------------------------
Edof=[1   1  2  3  4  5  6
      2   4  5  6  7  8  9
      3   7  8  9  10  11  12
      4   10  11  12  13  14  15
      5   13  14  15  16  17  18
      6   16  17  18  19  20  21
      7   19  20  21  22  23  24
      8   22  23  24  25  26  27
      9   25  26  27  28  29  30
      10   28  29  30  31  32  33
      11   31  32  33  34  35  36
      12   34  35  36  37  38  39
      13   37  38  39  40  41  42
      14   40  41  42  43  44  45
      15   43  44  45  46  47  48
      16   46  47  48  49  50  51
      17   49  50  51  52  53  54
      18   52  53  54  55  56  57
      19   55  56  57  58  59  60
      20   55  56  57  61  62  63
      21   61  62  63  109  110  111
      22   7  8  9  64  65  66
      23   64  65  66  67  68  69
      24   67  68  69  70  71  72
      25   70  71  72  76  77  78
      26   16  17  18  73  74  75
      27   73  74  75  76  77  78
      28   76  77  78  79  80  81
      29   79  80  81  85  86  87
      30   28  29  30  82  83  84
      31   82  83  84  85  86  87
      32   85  86  87  88  89  90
      33   88  89  90  94  95  96
      34   34  35  36  91  92  93
      35   91  92  93  94  95  96
      36   94  95  96  97  98  99
      37   97  98  99  103  104  105
      38   46  47  48  100  101  102
      39   100  101  102  103  104  105
      40   103  104  105  106  107  108
      41   106  107  108  109  110  111];
  
  lH=0.5475;
  err=1.e+6;
  i1=0;
  progression=0;
for H=lH
    
    % ------ list of coordinates  -----------------------------------
    Coord=[0  0; 0.1  0; 0.168  0; 0.2  0; 0.473  0; 0.778  0; 1.075  0; 1.175  0; 1.275  0; 1.388  0; 1.693  0; 1.998  0; 2.064  0; 2.164  0; 2.264  0; 2.608  0; 2.913  0; 3.093  0; 3.218  0; 3.296  0; 3.218  0.5*H; 0.168  0.5*H; 0.168  H; 0.473  H; 0.778  0.5*H; 0.778  H; 1.083  H; 1.388  0.5*H; 1.388  H; 1.693  H; 1.998  0.5*H; 1.998  H; 2.303  H; 2.608  0.5*H; 2.608  H; 2.913  H; 3.218  H];
    % ------ list of degrees-of-fredom  -----------------------------
    Dof=[1  2  3; 4  5  6; 7  8  9; 10  11  12; 13  14  15; 16  17  18; 19  20  21; 22  23  24; 25  26  27; 28  29  30; 31  32  33; 34  35  36; 37  38  39; 40  41  42; 43  44  45; 46  47  48; 49  50  51; 52  53  54; 55  56  57; 58  59  60; 61  62  63; 64  65  66; 67  68  69; 70  71  72; 73  74  75; 76  77  78; 79  80  81; 82  83  84; 85  86  87; 88  89  90; 91  92  93; 94  95  96; 97  98  99; 100  101  102; 103  104  105; 106  107  108; 109  110  111];
    % ------ generate element matrices, assemble in global matrices - 

   %% fixity factors
    lp2=0.1:0.01:1;
     p1=1;
    %f1obj=78.972; f2obj= 301.896; f3obj=318.916; f4obj=355.238; f5obj= 401.247; %5mm de largeur et 10mm de profondeur travée 1
    %f1obj=78.949; f2obj= 271.435; f3obj=312.662; f4obj=349.577; f5obj=398.422; %5mm de largeur et 20mm de profondeur travée 1
    %f1obj=79.011; f2obj= 301.227; f3obj=321.314; f4obj=357.796; f5obj= 398.357; %5mm de largeur et 10mm de profondeur travée 2
    f1obj=79.025; f2obj= 276.014; f3obj=317.696; f4obj=356.339; f5obj= 385.382; %5mm de largeur et 20mm de profondeur travée 2

    %%
     for p2=lp2
                    epv=[Ev Av Iv rho*Av 1.0 1.0]; %Barre de 20 à 35 et 38 à 39
                    epv36=[Ev Av Iv rho*Av 1.0 p1]; %Barre 36
                    epv37=[Ev Av Iv rho*Av p1 1.0]; %Barre 37
                    epv40=[Ev Av Iv rho*Av 1.0 p2]; %Barre 40
                    epv41=[Ev Av Iv rho*Av p2 1.0]; %Barre 41
                    eph=[Eh Ah Ih rho*Ah 1.0 1.0];  %Barre 1 à 19
                
        % ------ generate element matrices, assemble in global matrices - 
                K=zeros(111);     M=zeros(111);
                [Ex,Ey]=coordxtr(Edof,Coord,Dof,2);
                
            for i=20:35
              [k,m,c]=beam2d_sr(Ex(i,:),Ey(i,:),epv);
              K=assem(Edof(i,:),K,k);  M=assem(Edof(i,:),M,m);
            end
            
            for i=38:39
              [k,m,c]=beam2d_sr(Ex(i,:),Ey(i,:),epv);
              K=assem(Edof(i,:),K,k);  M=assem(Edof(i,:),M,m);
            end

            for i=1:19
              [k,m,c]=beam2d_sr(Ex(i,:),Ey(i,:),eph);
              K=assem(Edof(i,:),K,k);  M=assem(Edof(i,:),M,m);  
            end
            
            [k,m,c]=beam2d_sr(Ex(40,:),Ey(40,:),epv40);
            K=assem(Edof(40,:),K,k);  M=assem(Edof(40,:),M,m);
             
             [k,m,c]=beam2d_sr(Ex(41,:),Ey(41,:),epv41);
             K=assem(Edof(41,:),K,k);  M=assem(Edof(41,:),M,m);
              
              [k,m,c]=beam2d_sr(Ex(36,:),Ey(36,:),epv36);
              K=assem(Edof(36,:),K,k);  M=assem(Edof(36,:),M,m);
               
              [k,m,c]=beam2d_sr(Ex(37,:),Ey(37,:),epv37);
              K=assem(Edof(37,:),K,k);  M=assem(Edof(37,:),M,m);
            
    bc=[1 2 3 4 5 6 7 8 9 10 11 12 19 20 21 22 23 24 25 26 27 37 38 39 40 41 42 43 44 45 52 53 54 55 56 57 58 59 60]'; %apuis
        % ----- Eigenvalue analysis --------------------------------------

                [La,Egv]=eigen(K,M,bc);
                Freq1=sqrt(La(1))/(2*pi);
                Freq2=sqrt(La(2))/(2*pi);
                Freq3=sqrt(La(3))/(2*pi);
                Freq4=sqrt(La(4))/(2*pi);
                Freq5=sqrt(La(5))/(2*pi);
                Freq6=sqrt(La(6))/(2*pi);
                Freq7=sqrt(La(7))/(2*pi);
                Freq8=sqrt(La(8))/(2*pi);
                Freq9=sqrt(La(9))/(2*pi);
                    
               if imag(Freq1)==0
                  err1=sqrt((2*((f1obj-Freq1)^2)+2*((f2obj-Freq2)^2)+2*((f3obj-Freq3)^2)+2*((f4obj-Freq4)^2)+2*((f5obj-Freq5)^2))/10);
                  Freq_d_1=Freq1; Freq_d_2=Freq2 ; Freq_d_3=Freq3; Freq_d_4=Freq4; Freq_d_5=Freq5;
               elseif imag(Freq2)==0
                  err1=sqrt((2*((f1obj-Freq2)^2)+2*((f2obj-Freq3)^2)+2*((f3obj-Freq4)^2)+2*((f4obj-Freq5)^2)+2*((f5obj-Freq6)^2))/10);
                  Freq_d_1=Freq2; Freq_d_2=Freq3 ; Freq_d_3=Freq4; Freq_d_4=Freq5; Freq_d_5=Freq6;
               elseif imag(Freq3)==0
                  err1=sqrt((2*((f1obj-Freq3)^2)+2*((f2obj-Freq4)^2)+2*((f3obj-Freq5)^2)+2*((f4obj-Freq6)^2)+2*((f5obj-Freq7)^2))/10);
                  Freq_d_1=Freq3; Freq_d_2=Freq4 ; Freq_d_3=Freq5; Freq_d_4=Freq6; Freq_d_5=Freq7;
               elseif imag(Freq4)==0
                  err1=sqrt((2*((f1obj-Freq4)^2)+2*((f2obj-Freq5)^2)+2*((f3obj-Freq6)^2)+2*((f4obj-Freq5)^2)+2*((f5obj-Freq8)^2))/10);
                  Freq_d_1=Freq4; Freq_d_2=Freq5 ; Freq_d_3=Freq6; Freq_d_4=Freq7; Freq_d_5=Freq8;
               else
                  err1=sqrt((2*((f1obj-Freq5)^2)+2*((f2obj-Freq6)^2)+2*((f3obj-Freq7)^2)+2*((f4obj-Freq8)^2)+2*((f5obj-Freq9)^2))/10);
                  Freq_d_1=Freq5; Freq_d_2=Freq6 ; Freq_d_3=Freq7; Freq_d_4=Freq8; Freq_d_5=Freq9;
               end
                
                if (err1<err)
                    f1opt = Freq_d_1; f2opt=Freq_d_2 ;   f3opt=Freq_d_3; f4opt=Freq_d_4; f5opt=Freq_d_5;
                    p2opt=p2;
                    Hopt=H;
                    err=err1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                end
                i1=i1+1;
                progression=100*(i1/(length(lH)*length(lp1)))
     end
end

Hopt
p2opt
f1opt
f2opt
f3opt
f4opt
f5opt
err