function [Ke,Me,Ce]=beam3d_sr(ex,ey,ez,eo,ep)
% [Ke,Me]=beam2d_sr(ex,ey,ep)
% [Ke,Me,Ce]=beam2d_sr(ex,ey,ep)
%-------------------------------------------------------------
%   PURPOSE
%    Calculate the stiffness matrix Ke, the mass matrix Me 
%    and the damping matrix Ce for a 2D elastic Bernoulli
%    beam element with semi-rigid connections
%
%   INPUT:   ex = [x1 x2 x3]
%            ey = [y1 y2 y3]         element node coordinates 
%            eo = [xz yz zz];       orientation of local z axis
%
%            ep = [E G A Iy Iz Kv p1y p2y p1z p2z (a b)]
%                                 E: Young's modulus
%                                 G: Shear Modulus
%                                 A: cross section area
%                                 Iy : moment of inertia with respect to y
%                                 axis
%                                 Iz : moment of inertia with respect to z
%                                 axis                              
%                                 Kv: St Venant torsional stiffness
%                                 (Saint-Venant's theorem states that the simply connected cross section with maximal torsional rigidity is a circle. )
%								p1, p2 : fixity factor for origin and end
%                               a,b: damping coefficients,
%                                    Ce=aMe+bKe 
% 
%   OUTPUT:  Ke : element stiffness matrix (12 x 12)
%            Me : element mass matrix 
%            Ce : element damping matrix, optional
%-------------------------------------------------------------
% LAST MODIFIED: AM 5/05/2022
%-------------------------------------------------------------
    dx=ex(2)-ex(1);
    dy=ey(2)-ey(1);
    dz=ez(2)-ez(1);
    L=sqrt(dx^2+dy^2+dz^2);
    n1=[dx dy dz]/L;
    lc=sqrt(eo*eo'); 
    n3=eo/lc;
   

   %qX=0; qY=0; qZ=0; qW=0; 
   %if nargin>5; qX=eq(1); qY=eq(2); qZ=eq(3); qW=eq(4); end
   
   E=ep(1);    G=ep(2);    A=ep(3);     Iy=ep(4);  Iz=ep(5);  m=ep(6);  Kv = ep(7);
   p1y = ep(8); p2y = ep(9); p1z = ep(10); p2z = ep(11) ;
   a=0 ; b=0 ;
   if length(ep)==12 ; a=ep(11) ; b=ep(12) ; end
%
% b = Stiffness index
   b11 = 3*p1y/(4-(p1y*p2y));
   b12 = 3*p1y*p2y/(4-(p1y*p2y));
   b22 = 3*p2y/(4-(p1y*p2y));

   a11 = 3*p1z/(4-(p1z*p2z));
   a12 = 3*p1z*p2z/(4-(p1z*p2z));
   a22 = 3*p2z/(4-(p1z*p2z));

   Kle=[E*A/L   0                         0                        0            0                         0                          -E*A/L     0                           0                         0         0                       0;       
          0     4*E*Iz*(a11+a12+a22)/L^3  0                        0            0                         2*E*Iz*(a12+2*a11)/L^2       0       -4*E*Iz*(a11+a12+a22)/L^3    0                         0         0                       2*E*Iz*(a12+2*a22)/L^2;
          0     0                         4*E*Iy*(b11+b12+b22)/L^3 0           -2*E*Iy*(2*b11+b12)/L^2    0                            0        0                          -4*E*Iy*(b11+b12+b22)/L^3  0        -2*E*Iy*(2*b22+b12)/L^2  0; 
          0     0                         0                       G*Kv/L        0                         0                            0        0                           0                        -G*Kv/L    0                       0;
          0     0                        -2*E*Iy*(2*b11+b12)/L^2   0            4*E*Iy*b11/L              0                            0        0                          2*E*Iy*(2*b11+b12)/L^2     0         2*E*Iy*b12/L            0;
          0     2*E*Iz*(a12+2*a11)/L^2    0                        0            0                         4*E*Iz*a11/L                 0       -2*E*Iz*(a12+2*a11)/L^2      0                         0         0                       2*E*Iz*a12/L;
         -E*A/L 0                         0                        0            0                         0                           E*A/L     0                           0                         0         0                       0;
          0    -4*E*Iz*(a11+a12+a22)/L^3  0                        0            0                        -2*E*Iz*(a12+2*a11)/L^2       0        4*E*Iz*(a11+a12+a22)/L^3    0                         0         0                      -2*E*Iz*(a12+2*a22)/L^2;
          0     0                        -4*E*Iy*(b11+b12+b22)/L^3 0            2*E*Iy*(2*b11+b12)/L^2    0                            0        0                          4*E*Iy*(b11+b12+b22)/L^3   0         2*E*Iy*(2*b22+b12)/L^2  0;
          0     0                         0                      -G*Kv/L        0                         0                            0        0                           0                         G*Kv/L    0                       0;
          0     0                        -2*E*Iy*(2*b22+b12)/L^2   0            2*E*Iy*b12/L              0                            0        0                          2*E*Iy*(2*b22+b12)/L^2     0         4*E*Iy*b22/L            0;
          0     2*E*Iz*(a12+2*a22)/L^2    0                        0            0                         2*E*Iz*a12/L                 0       -2*E*Iz*(a12+2*a22)/L^2      0                         0         0                       4*E*Iz*a22/L];
%


   dz = 4-p1z*p2z;
   f1z=32*p1z^2*p2z^2-55*p1z^2*p2z+32*p1z^2+50*p1z*p2z^2-328*p1z*p2z+224*p1z+32*p2z^2-196*p2z+560;
   f2z=25*p1z^2*p2z^2-86*p1z^2*p2z+64*p1z^2+32*p1z*p2z^2-160*p1z*p2z+224*p1z;
   f3z=41*p1z^2*p2z^2+5*p1z^2*p2z-64*p1z^2+5*p1z*p2z^2-184*p1z*p2z-28*p1z-64*p2z^2-28*p2z+560;
   f4z=55*p1z^2*p2z^2-64*p1z^2*p2z-38*p1z*p2z^2-100*p1z*p2z-128*p2z^2+392*p2z;
   f5z=32*p1z^2-31*p1z^2*p2z+8*p1z^2*p2z^2;
   f6z=31*p1z^2*p2z^2-64*p1z^2*p2z-64*p1z*p2z^2+124*p1z*p2z;

   dy = 4-p1y*p2y;
   f1y= 32*p1y^2*p2y^2-55*p1y^2*p2y+32*p1y^2+50*p1y*p2y^2-328*p1y*p2y+224*p1y+32*p2y^2-196*p2y+560;
   f2y= 25*p1y^2*p2y^2-86*p1y^2*p2y+64*p1y^2+32*p1y*p2y^2-160*p1y*p2y+224*p1y;
   f3y= 41*p1y^2*p2y^2+5*p1y^2*p2y-64*p1y^2+5*p1y*p2y^2-184*p1y*p2y-28*p1y-64*p2y^2-28*p2y+560;
   f4y= 55*p1y^2*p2y^2-64*p1y^2*p2y-38*p1y*p2y^2-100*p1y*p2y-128*p2y^2+392*p2y;
   f5y= 32*p1y^2-31*p1y^2*p2y+8*p1y^2*p2y^2;
   f6y= 31*p1y^2*p2y^2-64*p1y^2*p2y-64*p1y*p2y^2+124*p1y*p2y;


   Mle= (m*L/420)*[
				140       0         0       0        0           0  	   70        0        0        0         0           0   ;
                 0   	4*f1z/dz^2  0       0        0     2*L*f2z/dz^2     0    2*f3z/dz^2   0        0         0     -L*f4z/dz^2 ;
                 0        0     4*f1y/dy^2  0  -2*L*f2y/dy^2     0          0        0   2*f3y/dy^2    0    L*f4y/dy^2       0   ;
                 0        0         0    140*(Iz/A)  0           0          0        0        0     -70*(Iz/A)   0           0   ;
                 0        0  -2*L*f2y/dy^2  0  4*L^2*f5y/dy^2    0          0        0   -L*f4y/dy^2   0   -L^2*f6y/dy^2     0   ;     
                 0    2*L*f2z/dz^2  0       0        0   4*L^2*f5z/dz^2     0   L*f4z/dz^2    0        0         0      -L^2*f6z/dz^2;
                70        0         0       0        0           0  	   140       0        0        0         0           0   ;
                 0    2*f3z/dz^2    0       0        0     L*f4z/dz^2   	0   4*f1z/dz^2    0        0         0     -2*L*f2z/dz^2 ;
                 0   	  0     2*f3y/dy^2  0   -L*f4y/dy^2      0   	    0        0   4*f1y/dy^2    0    2*L*f2y/dy^2     0   ;
                 0        0         0   -70*(Iz/A)   0           0          0        0        0     140*(Iz/A)   0           0   ;
                 0        0     L*f4y/dy^2  0   -L^2*f6y/dy^2    0          0        0  2*L*f2y/dy^2   0   4*L^2*f5y/dy^2    0   ;
                 0    -L*f4z/dz^2   0       0        0    -L^2*f6z/dz^2     0  -2*L*f2z/dz^2  0        0         0     4*L^2*f5z/dz^2;
                 
				 ];
%

    Cle=a*Mle+b*Kle;

  n2(1)= n3(2)*n1(3)-n3(3)*n1(2);
  n2(2)= n3(3)*n1(1)-n3(1)*n1(3);
  n2(3)= n3(1)*n1(2)-n1(1)*n3(2);
 
%
  
%
  G=[ n1(1)     n1(2)   n1(3)  0      0       0     0      0       0      0      0      0;
      n2(1)     n2(2)   n2(3)  0      0       0     0      0       0      0      0      0;
      n3(1)     n3(2)   n3(3)  0      0       0     0      0       0      0      0      0;
       0         0        0   n1(1)  n1(2)   n1(3)  0      0       0      0      0      0;
       0         0        0   n2(1)  n2(2)   n2(3)  0      0       0      0      0      0;
       0         0        0   n3(1)  n3(2)   n3(3)  0      0       0      0      0      0;
       0         0        0   0       0       0   n1(1)   n1(2)   n1(3)   0      0      0;
       0         0        0   0       0       0   n2(1)   n2(2)   n2(3)   0      0      0;
       0         0        0   0       0       0   n3(1)   n3(2)   n3(3)   0      0      0;
       0         0        0   0       0       0     0      0       0    n1(1)   n1(2)  n1(3);    
       0         0        0   0       0       0     0      0       0    n2(1)   n2(2)  n2(3); 
       0         0        0   0       0       0     0      0       0    n3(1)   n3(2)  n3(3)];
%
    Ke=G'*Kle*G;  Me=G'*Mle*G; Ce=G'*Cle*G;
%--------------------------end--------------------------------

 
