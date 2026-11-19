// === Paramètres ===
L = 5.0;        // Longueur de la demi-éprouvette
W = 2.0;        // Largeur de la demi-éprouvette
LL = 2.0;
WW = 1.0;
T = 0.5;
hh = 0.25;        // Taille cible des éléments
h = .5;        // Taille cible des éléments
R = 2.;
Delta = W-WW;
Theta = Acos(Delta/R);
LLL = LL+R*Sin(Theta);

// === Points ===
Point(1) = {0, 0, 0, hh};     // bas gauche
Point(2) = {L, 0, 0, hh};     // bas droit
Point(3) = {L, WW, 0, hh};     // haut droit
Point(4) = {L-LL, WW, 0, hh};
Point(5) = {L-LLL, W, 0, hh};     
Point(6) = {0, W, 0, h};     // haut gauche
Point(7) = {L-LL, WW+R, 0, hh}; // centre rayon
Point(8) = {0, 0, T, hh};     // bas gauche
Point(9) = {L, 0, T, hh};     // bas droit
Point(10) = {L, WW, T, hh};     // haut droit
Point(11) = {L-LL, WW, T, hh};
Point(12) = {L-LLL, W, T, hh};     
Point(13) = {0, W, T, h};     // haut gauche
Point(14) = {L-LL, WW+R, T, hh}; // centre rayon


// === Lignes ===
Line(1) = {1, 2};  // bas
Line(2) = {2, 3};  // droite
Line(3) = {3, 4};  // haut
Line(4) = {5, 6};  // haut
Line(5) = {6, 1};  // gauche
Circle(6) = {4, 7, 5};
Line(7) = {8, 9};  // bas
Line(8) = {9, 10};  // droite
Line(9) = {10, 11};  // haut
Line(10) = {12, 13};  // haut
Line(11) = {13, 8};  // gauche
Circle(12) = {11, 14, 12};
Line(13) = {1, 8};
Line(14) = {2, 9}; 
Line(15) = {3, 10};
Line(16) = {4, 11};
Line(17) = {5, 12};
Line(18) = {6, 13};
// === Contours ===
Curve Loop(1) = {1, 2, 3, 4, 5, 6}; // face inf
Curve Loop(2) = {7, 8, 9, 10, 11, 12}; // face sup
Curve Loop(3) = {18, 11, -13, -5}; // face inf
Curve Loop(4) = {15, -8, -14, 2}; // face sup
Curve Loop(5) = {17, 10, -18, -4};
Curve Loop(6) = {15, 9, -16, -3};
Curve Loop(7) = {14, -7, -13, 1};
Curve Loop(8) = {16, 12, -17, -6};
// === Surface ===
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Surface(8) = {8};

Surface Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Volume(1) = {1};

Physical Surface("ZSYM") = {1};
Physical Surface("XDSP") = {3};
Physical Surface("YSYM") = {7};
Physical Surface("XSYM") = {4};



Physical Volume("Specimen") = {1};
// --- Options de maillage ---
Mesh.ElementOrder = 2;
Mesh.RecombineAll = 0;
Mesh.SaveAll = 0; 
Mesh 3;

