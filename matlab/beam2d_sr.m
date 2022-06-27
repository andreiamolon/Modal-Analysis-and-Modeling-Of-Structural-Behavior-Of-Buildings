function [Ke,Me,Ce]=beam2d_sr(ex,ey,ep)
% [Ke,Me]=beam2d_sr(ex,ey,ep)
% [Ke,Me,Ce]=beam2d_sr(ex,ey,ep)
%-------------------------------------------------------------
%   PURPOSE
%    Calculate the stiffness matrix Ke, the mass matrix Me 
%    and the damping matrix Ce for a 2D elastic Bernoulli
%    beam element with semi-rigid connections
%
%   INPUT:   ex = [x1 x2]
%            ey = [y1 y2]         element node coordinates 
%
%            ep = [E A I m p1 p2(a b)]  
%                                 E: Young's modulus
%                                 A: cross section area
%                                 I: moment of inertia
%                                 m: mass per unit length
%								p1, p2 : fixity factor for origin and end
%                               a,b: damping coefficients,
%                                    Ce=aMe+bKe 
% 
%   OUTPUT:  Ke : element stiffness matrix (6 x 6)
%            Me : element mass matrix 
%            Ce : element damping matrix, optional
%-------------------------------------------------------------
% LAST MODIFIED: MM 18/02/2013
%-------------------------------------------------------------
   b=[ ex(2)-ex(1); ey(2)-ey(1) ];
   L=sqrt(b'*b);  n=b/L;
   
   E=ep(1);    A=ep(2);    I=ep(3);     m=ep(4); p1=ep(5);  p2 = ep(6);
   a=0 ; b=0 ;
   if length(ep)==8 ; a=ep(7) ; b=ep(8) ; end
%
   b11 = 3*p1/(4-p1*p2);
   b12 = 3*p1*p2/(4-p1*p2);
   b22 = 3*p2/(4-p1*p2);
   Kle=[E*A/L   0                         0                      -E*A/L      0                         0 ;
          0     4*E*I*(b11+b12+b22)/L^3   2*E*I*(2*b11+b12)/L^2  0           -4*E*I*(b11+b12+b22)/L^3  2*E*I*(b12+2*b22)/L^2;
          0     2*E*I*(2*b11+b12)/L^2     4*E*I*b11/L            0           -2*E*I*(2*b11+b12)/L^2     2*E*I*b12/L;
        -E*A/L  0                         0                      E*A/L       0                         0 ;
          0     -4*E*I*(b11+b12+b22)/L^3  -2*E*I*(2*b11+b12)/L^2  0           4*E*I*(b11+b12+b22)/L^3   -2*E*I*(b12+2*b22)/L^2;
          0     2*E*I*(b12+2*b22)/L^2     2*E*I*b12/L            0           -2*E*I*(b12+2*b22)/L^2    4*E*I*b22/L];
%
   d=4-p1*p2;
   f1=560+224*p1+32*p1*p1-196*p2-328*p1*p2-55*p1*p1*p2+32*p2*p2+50*p1*p2*p2+32*p1*p1*p2*p2;
   f2=224*p1+64*p1*p1-160*p1*p2-86*p1*p1*p2+32*p1*p2*p2+25*p1*p1*p2*p2;
   f3=560-28*p1-64*p1*p1-28*p2-184*p1*p2+5*p1*p1*p2-64*p2*p2+5*p1*p2*p2+41*p1*p1*p2*p2;
   f4=392*p2-100*p1*p2-64*p1*p1*p2-128*p2*p2-38*p1*p2*p2+55*p1*p1*p2*p2;
   f5=32*p1*p1-31*p1*p1*p2+8*p1*p1*p2*p2;
   f6=124*p1*p2-64*p1*p1*p2-64*p1*p2*p2+31*p1*p1*p2*p2;
   
   Mle=m*L/420/d/d*[
				140*d*d   0     0    	  70*d*d    0      0    ;
                 0   	4*f1    2*L*f2    0   		2*f3  -L*f4 ;
                 0   	2*L*f2  4*L^2*f5  0  		L*f4  -L^2*f6 ;
                70*d*d    0     0   	  140*d*d   0      0   ;
                 0   	2*f3    L*f4   	  0         4*f1  -2*L*f2 ;
                 0      -L*f4   -L^2*f6   0        -2*L*f2 4*L^2*f5
				 ];
%
   Cle=a*Mle+b*Kle;
%
   G=[n(1) n(2)  0    0    0   0;
     -n(2) n(1)  0    0    0   0;
       0    0    1    0    0   0;
       0    0    0   n(1) n(2) 0;
       0    0    0  -n(2) n(1) 0;
       0    0    0    0    0   1];
%
   Ke=G'*Kle*G;    Me=G'*Mle*G;    Ce=G'*Cle*G;
%--------------------------end--------------------------------
 
