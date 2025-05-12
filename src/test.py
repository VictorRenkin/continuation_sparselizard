import numpy as np
import csv

def lire_array_depuis_fichier(nom_fichier):
    """Lit un fichier .txt et retourne un array numpy à partir d'une ligne de valeurs séparées par des virgules."""
    try:
        with open(nom_fichier, 'r') as f:
            ligne = f.readline()
            valeurs = [float(x.strip()) for x in ligne.split(',')]
            return np.array(valeurs)
    except Exception as e:
        print(f"Erreur lors de la lecture du fichier : {e}")
        return None

array_x_2 = lire_array_depuis_fichier("X_2.txt")
array_y_2 = lire_array_depuis_fichier("B.txt")

array_x_2 = array_x_2[1:]
array_y_2 = array_y_2[1:]

print('length array_x_2', len(array_x_2))
print('length array_y_2', len(array_y_2))

if len(array_x_2) == len(array_y_2):
    # Écriture dans un fichier CSV
    with open("comparaison.csv", "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["X_2", "B"])  # En-têtes de colonne
        for x, y in zip(array_x_2, array_y_2):
            writer.writerow([x, y])
    print("Fichier CSV 'comparaison.csv' généré avec succès.")
    
    produit_scalaire = np.dot(array_x_2, array_y_2)
    print("Produit scalaire :", produit_scalaire)
else:
    print("Erreur : Les deux tableaux n'ont pas la même longueur.")
