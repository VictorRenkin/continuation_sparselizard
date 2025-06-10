// Fichier : ligne_2_elements_force.geo

// Trois points alignés
Point(1) = {0, 0, 0, 1.0};   // Clamped
Point(2) = {1, 0, 0, 1.0};   // Force
Point(3) = {2, 0, 0, 1.0};   // Clamped

// Deux lignes
Line(1) = {1, 2};
Line(2) = {2, 3};

// Forcer le nombre de divisions : 2 éléments => 3 nœuds
Transfinite Line{1} = 2; // 2 points => 1 élément
Transfinite Line{2} = 2; // 2 points => 1 élément

// Indiquer que les points sont partie d’un maillage transfinite
Transfinite Point{1, 2, 3};

// Physicals
Physical Point(2) = {1, 3};
Physical Point(3) = {2};
Physical Line(1) = {1, 2};
