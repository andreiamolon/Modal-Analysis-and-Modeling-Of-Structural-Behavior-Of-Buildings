
%----------------------------------------------------------------
% PURPOSE 
%    Set up the fe-model and perform eigenvalue analysis
%    for a simple frame structure.
%----------------------------------------------------------------

% REFERENCES
%     G"oran Sandberg 1994-03-08 
%     Karl-Gunnar Olsson 1995-09-29
%----------------------------------------------------------------
 echo on

% ------ Generate the model ------------------------------------------

% ------ material data ------------------------------------------
E=3e10;                 rho=2500;
Av=0.36e-2;           Iv=1.08e-6;             % Poteaux et poutres
Ah=0.0375;           Ih=1.953125e-4;            % Longrines
epv=[E Av Iv rho*Av];   eph=[E Ah Ih rho*Ah];
H=0.6004;
p=0.5;      %Facteur de fixité

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
      40   103  104  105  112  113  114
      41   115  116  117  109  110  111
      42   112  113  114  106  107  108
      43   106  107  108  115  116  117];
% ------ list of coordinates  -----------------------------------
Coord=[0  0; 0.1  0; 0.168  0; 0.2  0; 0.473  0; 0.778  0; 1.075  0; 1.175  0; 1.275  0; 1.388  0; 1.693  0; 1.998  0; 2.064  0; 2.164  0; 2.264  0; 2.608  0; 2.913  0; 3.093  0; 3.218  0; 3.296  0; 3.218  0.5*H; 0.168  0.5*H; 0.168  H; 0.473  H; 0.778  0.5*H; 0.778  H; 1.083  H; 1.388  0.5*H; 1.388  H; 1.693  H; 1.998  0.5*H; 1.998  H; 2.303  H; 2.608  0.5*H; 2.608  H; 2.913  H; 3.218  H; 2.9105 H; 2.9155 H];
% ------ list of degrees-of-fredom  -----------------------------
Dof=[1  2  3; 4  5  6; 7  8  9; 10  11  12; 13  14  15; 16  17  18; 19  20  21; 22  23  24; 25  26  27; 28  29  30; 31  32  33; 34  35  36; 37  38  39; 40  41  42; 43  44  45; 46  47  48; 49  50  51; 52  53  54; 55  56  57; 58  59  60; 61  62  63; 64  65  66; 67  68  69; 70  71  72; 73  74  75; 76  77  78; 79  80  81; 82  83  84; 85  86  87; 88  89  90; 91  92  93; 94  95  96; 97  98  99; 100  101  102; 103  104  105; 106  107  108; 109  110  111; 112 113 114; 115 116 117];
% ------ generate element matrices, assemble in global matrices - 
K=zeros(117);     M=zeros(117);
[Ex,Ey]=coordxtr(Edof,Coord,Dof,2);
for i=20:41
  [k,m,c]=beam2d(Ex(i,:),Ey(i,:),epv);
  K=assem(Edof(i,:),K,k);  M=assem(Edof(i,:),M,m);  
end
for i=1:19
  [k,m,c]=beam2d(Ex(i,:),Ey(i,:),eph);
  K=assem(Edof(i,:),K,k);  M=assem(Edof(i,:),M,m);  
end
for i=42:43
  [k,m,c]=beam2d(Ex(i,:),Ey(i,:),epv);
  K=assem(Edof(i,:),K,p*k);  M=assem(Edof(i,:),M,p*m);  
end
 

% ----- Draw a plot of the element mesh --------------------------


clf;     eldraw2(Ex,Ey,[1 2 2],Edof);
grid;    title('2-D Frame Structure') 

% ----- Eigenvalue analysis --------------------------------------

b=[1 2 3 4 5 6 7 8 9 10 11 12 19 20 21 22 23 24 25 26 27 37 38 39 40 41 42 43 44 45 52 53 54 55 56 57 58 59 60]';
[La,Egv]=eigen(K,M,b);
Freq=sqrt(La)/(2*pi);


% ----- plot one eigenmode ---------------------------------------

figure(1),    clf,     grid,     title('The second eigenmode'), 
eldraw2(Ex,Ey,[2 3 1]); 
Edb=extract(Edof,Egv(:,2));      eldisp2(Ex,Ey,Edb,[1 2 2]);
FreqText=num2str(Freq(1));       text(.5,1.75,FreqText);

% ----- plot one eigenmode ---------------------------------------

figure(4),    clf,     grid,     title('The first eigenmode'), 
eldraw2(Ex,Ey,[2 3 1]); 
Edb=extract(Edof,Egv(:,1));      eldisp2(Ex,Ey,Edb,[1 2 2]);
FreqText=num2str(Freq(1));       text(.5,1.75,FreqText);


% ----- plot eight eigenmodes ------------------------------------

figure(2), clf, axis('equal'), hold on, axis off
sfac=0.5;
title('The first eigenmode (Hz)' )
for i=1:4
  Edb=extract(Edof,Egv(:,i));
  Ext=Ex+(i-1)*3;                eldraw2(Ext,Ey,[2 3 1]); 
  eldisp2(Ext,Ey,Edb,[1 2 2],sfac);
  FreqText=num2str(Freq(i));     text(3*(i-1)+.5,1.5,FreqText);
end
%Eyt=Ey-4; 
%for i=5:8
  %Edb=extract(Edof,Egv(:,i));
  %Ext=Ex+(i-5)*3;                eldraw2(Ext,Eyt,[2 3 1]); 
  %eldisp2(Ext,Eyt,Edb,[1 2 2],sfac);
  %FreqText=num2str(Freq(i));     text(3*(i-5)+.5,-2.5,FreqText);
%end

% -------------------- end ---------------------------------------
echo off
