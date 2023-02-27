syms q2 DG DN DC dG dN dS kplus kminus N0 G0 H_GS H_NS

LinearizedRDMatrix = -[(dG+kplus*N0+DG*q2) (kplus*N0) (-kplus*N0) (-dG*H_GS) ;
    (kplus*G0) (dN+kplus*G0+DN*q2) (-kplus*G0) (-dN*H_NS);
    (-kminus) (-kminus) (kminus + DC*q2) (0);
    (-dS) (0) (0) (dS)];

determinant_value = det(-LinearizedRDMatrix);

polynomial_coefficients = fliplr(coeffs(determinant_value,q2));

b3 = simplify(polynomial_coefficients(1)) % b3 >0 for all params
b2 = simplify(polynomial_coefficients(2)) % b2 > 0 given H_GS < 0
b1 = simplify(polynomial_coefficients(3)) % b1 < 0 only if H_NS < 0
b0 = simplify(polynomial_coefficients(4)) % b0 > 0 given H_GS < 0

%% Include degradation of complex

syms q2 DG DN DC dG dN dS kplus N0 G0 H_GS H_NS f keff

LinearizedRDMatrix = -[(dG+kplus*N0+DG*q2) (kplus*N0) (-kplus*N0*(1-f)) (-(dG+f*kplus*N0)*H_GS) ;
    (kplus*G0) (dN+kplus*G0+DN*q2) (-kplus*G0*(1-f)) (-(dN+f*kplus*G0)*H_NS);
    (-keff) (-keff) (keff + DC*q2) (0);
    (-dS) (0) (0) (dS)];

determinant_value = det(-LinearizedRDMatrix);

polynomial_coefficients = fliplr(coeffs(determinant_value,q2));

b3 = simplify(polynomial_coefficients(1)) % b3 >0 for all params
b2 = simplify(polynomial_coefficients(2)) % b2 > 0 given H_GS < 0
b1 = simplify(polynomial_coefficients(3)) % b1 < 0 only if H_NS < 0
b0 = simplify(polynomial_coefficients(4)) % b0 > 0 by det(-F^{RD}(q=0)) > 0 condition