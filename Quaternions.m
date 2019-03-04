function [Qtd,A,Eul] = Quaternions(PQR,Qt)
% Date: 03/04/2019
% Author: Leland Klein

% Qtd: the Quaternion derivatives from current Quaternions and body rates
    p 		= PQR(1);
	q 		= PQR(2);
	r 		= PQR(3);
	q1		= Qt(1);
	q2		= Qt(2);
	q3 		= Qt(3);
	q4 		= Qt(4); 
    
    q1dot = 0.5 * (  q4 * p - q3 * q + q2 * r );
	q2dot = 0.5 * (  q3 * p + q4 * q - q1 * r );
	q3dot = 0.5 * ( -q2 * p + q1 * q + q4 * r );
	q4dot = 0.5 * ( -q1 * p - q2 * q - q3 * r );
    
    Qtd =[ q1dot q2dot  q3dot q4dot];
  
% A: the cosine matrix from current Quaternions
    
    a11 = q1^2 - q2^2 - q3^2 + q4^2;
	a12 = 2 * ( q1 * q2 + q3 * q4 );
	a13 = 2 * ( q1 * q3 - q2 * q4 );
	a21 = 2 * ( q1 * q2 - q3 * q4 );
	a22 = q2^2 - q1^2 - q3^2 + q4^2;
	a23 = 2 * ( q2 * q3 + q1 * q4 );
	a31 = 2 * ( q1 * q3 + q2 * q4 );
	a32 = 2 * ( q2 * q3 - q1 * q4 );
	a33 = q3^2 - q1^2 - q2^2 + q4^2;

	A =    [a11 a12 a13;
	 		a21 a22 a23;
	 		a31 a32 a33];

% Eul: the Euler angles from current Quaternions
 	
	 	psi 	= atan2( a12, a11 );
		theta 	= atan2( -a13, sqrt( a11 * a11 + a12 * a12 ) );
		phi 	= atan2( a23, a33 );
        
    Eul =[phi  theta psi];
end

