import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# Définition globale des paramètres de police et de taille pour tous les graphiques
plt.rc('font', family='serif')  # Police avec empattements, comme Times
plt.rc('text', usetex=True)  # Utiliser LaTeX pour le texte dans les figures
plt.rcParams.update({
    'font.size': 18,         # Taille générale de la police
    'legend.fontsize': 23,   # Taille de la police des légendes
    'axes.labelsize': 20,    # Taille de la police des étiquettes des axes
    'axes.titlesize': 14,    # Taille de la police pour les titres des axes
    'xtick.labelsize': 16,   # Taille de la police pour les graduations de l'axe des X
    'ytick.labelsize': 16,   # Taille de la police pour les graduations de l'axe des Y
    'text.usetex': True,     # Confirmer l'utilisation de LaTeX pour tout le texte
    'figure.titlesize': 18   # Taille de la police pour les titres des figures
})

color_list = [
    "#6e343d", "#007070", "#f07f3c", "#5b57a2", "#7db928", "#e62d31",
    "#005ca9", "#00843b", "#f8aa00", "#5b257d", "#8c8b82"
]


def extract_value_and_time(df, value_column, time_column="Time"):
    """
    Extrait les paires (valeur, temps) valides à partir d'une colonne du DataFrame,
    où les valeurs sont numériques (i.e., différentes de '-').
    
    Retourne :
        - values: vecteur numpy des valeurs valides
        - times: vecteur numpy des temps associés
    """
    values = pd.to_numeric(df[value_column], errors='coerce')
    times = pd.to_numeric(df[time_column], errors='coerce')
    valid_mask = values.notna()

    return values[valid_mask].to_numpy(), times[valid_mask].to_numpy()

def smart_concatenate(base_array, new_array):
    """
    Concatène intelligemment new_array à base_array pour assurer une continuité naturelle :
    - Si new_array commence proche de la fin de base_array, on le garde tel quel.
    - Si new_array commence proche du début de base_array, on l'inverse avant de le concaténer.
    !!!!!!!!!!!!! A l'envers si c'est un backward !!!!!!!!!!!!!!
    """
    if len(base_array) == 0:
        return new_array

    # Calculer la distance entre les extrémités
    dist_start = abs(base_array[-1] - new_array[0])
    dist_end = abs(base_array[-1] - new_array[-1])

    if dist_start <= dist_end:
        return np.concatenate([base_array, new_array])
    else:
        return np.concatenate([base_array, new_array[::-1]])

df = pd.read_csv('output_Pseudo-arcegth rapport mesh.csv')
# Concaténation des vecteurs pour obtenir un seul vecteur frequency et forward
forward_1, frequency_1 = extract_value_and_time(df, "Forward-1")
forward_2, frequency_2 = extract_value_and_time(df, "Forward-2")
forward_3, frequency_3 = extract_value_and_time(df, "Forward-3")

predictor_forward_1, predictor_frequency_1 = extract_value_and_time(df, "Predictor-1")
predictor_forward_2, predictor_frequency_2 = extract_value_and_time(df, "Predictor-2")
predictor_forward_3, predictor_frequency_3 = extract_value_and_time(df, "Predictor-3")
# Application au vecteur fréquence
frequency = smart_concatenate(frequency_1, frequency_2)
frequency = smart_concatenate(frequency, frequency_3)

# # Et pour le vecteur forward correspondant
forward = smart_concatenate(forward_1, forward_2)
forward = smart_concatenate(forward, forward_3)
bifurcation, frequency_bif = extract_value_and_time(df, "Bifurcation")

predictor_frequency = smart_concatenate(predictor_frequency_1, predictor_frequency_2)
predictor_frequency = smart_concatenate(predictor_frequency, predictor_frequency_3)

predictor_forward = smart_concatenate(predictor_forward_1, predictor_forward_2)
predictor_forward = smart_concatenate(predictor_forward, predictor_forward_3)
print("Len predictor", len(predictor_forward))
print("Len normal", len(forward))




# # Données downward
# x_down = [
#     1025.8992805755395, 1026.2589928057555, 1026.6187050359713, 1028.0575539568345,
#     1029.136690647482, 1031.6546762589928, 1034.5323741007194, 1037.0503597122301,
#     1039.568345323741, 1042.8057553956835, 1046.7625899280574, 1049.6402877697842,
#     1053.2374100719423, 1056.8345323741007, 1060.7913669064749, 1065.4676258992806,
#     1070.1438848920864, 1075.5395683453237, 1080.9352517985612, 1087.7697841726617,
#     1094.9640287769785, 1102.158273381295, 1112.9496402877699, 1142.0863309352517,
#     1180.5755395683454, 1215.4676258992806, 1268.3453237410072
# ]
# x_down = [x / (2 * np.pi) for x in x_down]
# y_down = [
#     0.003658536585365854, 0.0033170731707317077, 0.002878048780487805,
#     0.0024878048780487805, 0.0021219512195121953, 0.001829268292682927,
#     0.0015853658536585367, 0.0014146341463414636, 0.0012195121951219512,
#     0.0010975609756097562, 0.000951219512195122, 0.0008048780487804879,
#     0.0007073170731707318, 0.0006341463414634147, 0.0005853658536585367,
#     0.0005121951219512195, 0.00043902439024390245, 0.00041463414634146346,
#     0.00036585365853658537, 0.0003414634146341464, 0.0002682926829268293,
#     0.0002682926829268293, 0.00024390243902439024, 0.00014634146341463417,
#     0.00004878048780487805, 0.00004878048780487805, 0
# ]
plt.figure(figsize=(12, 8))  # Augmenter la taille de la figure pour une meilleure visibilité
plt.scatter(frequency_bif, bifurcation, color=color_list[1], label=r'Bifurcation', marker='x', s=60, zorder=3)  # Augmenter la taille des marqueurs
for i in range(len(predictor_forward)):
    x = [frequency[i], predictor_frequency[i]]
    y = [forward[i], predictor_forward[i]]
    plt.plot(x, y, color=color_list[3], linestyle='-', marker='o',
         linewidth=1, markersize=4)
plt.plot(forward[i]+1000, predictor_forward[i]+10000, color=color_list[3], label=r'Predictor')
plt.plot(frequency, forward, color=color_list[0], linestyle='-', marker='o',
         linewidth=1, markersize=4 ,label=r'Forward', zorder=2)



plt.xlabel(r"Frequency [Hz]")
plt.ylabel(r"Max displacement [m]")
plt.xlim(np.min(frequency), 185)  # Ajouter une marge pour mieux voir les extrémités
plt.ylim(0, np.max(forward) * 1.2)  # Ajouter une marge pour mieux voir les données
plt.legend(loc='best')  # Ajuster la taille de la légende
plt.tight_layout()  # Optimiser l'agencement des éléments
plt.savefig("NLFR_enhanced.png", bbox_inches='tight', dpi=300)



df = pd.read_csv('cantilaveur.csv')
# # Concaténation des vecteurs pour obtenir un seul vecteur frequency et forward
backward_1, frequency_1 = extract_value_and_time(df, "Backward-1")
backward_2, frequency_2 = extract_value_and_time(df, "Backward-2")
backward_3, frequency_3 = extract_value_and_time(df, "Backward-3")

df_txt = pd.read_csv("node_0.005000_0.050000_1.000000.txt", sep=';')

df_txt["amp_norm"] = np.sqrt(df_txt["ampx"]**2 + df_txt["ampy"]**2 + df_txt["ampz"]**2)




backward = smart_concatenate(backward_3, backward_2)
backward = smart_concatenate(backward, backward_1)
frequency = smart_concatenate(frequency_3, frequency_2)
frequency = smart_concatenate(frequency, frequency_1)
bifurcation, frequency_bif = extract_value_and_time(df, "Bifurcation")
plt.figure(figsize=(12, 8))  # Augmenter la taille de la figure pour une meilleure visibilité
plt.plot(frequency, backward, color=color_list[0], marker='o', markersize=4, label=r'Backward', linewidth=1.5, zorder=2)  # Augmenter la taille des marqueurs et l'épaisseur de la ligne
plt.plot(df_txt["freq"]/2/np.pi, df_txt["amp_norm"],color = color_list[2], marker='o', markersize=4, label=r'Blahos ', linewidth=1.5, zorder=2)  # Augmenter la taille des marqueurs et l'épaisseur de la ligne
plt.scatter(frequency_bif, bifurcation, color=color_list[1], label=r'Bifurcation', marker='x', s=60, zorder=3)  # Augmenter la taille des marqueurs
plt.xlim(3.9,4.05)
plt.ylim(0.1, np.max(backward) * 1.1)  # Ajouter une marge pour mieux voir les données
plt.xlabel(r"Frequency [Hz]")
plt.ylabel(r"Max displacement [m]")
plt.legend(loc='best')  # Ajouter la légende
plt.tight_layout()  # Optimiser l'agencement des éléments
plt.savefig("NLFR_backward_cantivilateur.pdf", bbox_inches='tight', dpi=300)


df_MPI = pd.read_csv('MPI_Timing_Data.csv')


# Préparer les données
x_labels = df_MPI['MPI'].astype(str)
generate_times = df_MPI['time_generate']/60
solve_times = df_MPI['time_solve']/60
time_other = df_MPI['time_all']/60 - (generate_times + solve_times)



eff_general_times = generate_times/generate_times[0] * 100
eff_solve_times = solve_times/solve_times[0] * 100
eff_other_times = time_other/time_other[0] * 100

# Créer un graphique en bâtons empilés
plt.figure(figsize=(10, 6))
plt.bar(x_labels, generate_times, label=r'Jacobian setup', color=color_list[0], edgecolor='black')
plt.bar(x_labels, solve_times, bottom=generate_times, label=r'Solve', color=color_list[1], edgecolor='black')
plt.bar(x_labels, time_other, bottom=generate_times + solve_times, label=r'Other', color=color_list[2], edgecolor='black')

plt.xlabel(r"MPI [-]")
plt.ylabel(r"Time [min]")
plt.legend()
plt.tight_layout()
plt.savefig("MPI_Timing_Data.pdf", bbox_inches='tight', dpi=300)
plt.close()

plt.figure(figsize=(10, 6))
plt.plot(x_labels, eff_general_times, label=r'Jacobian setup', color=color_list[0], marker='o')
plt.plot(x_labels, eff_solve_times, label=r'Solve', color=color_list[1], marker='s')
# plt.plot(x_labels, eff_other_times, label=r'Other', color=color_list[2], marker='^')
plt.xlabel(r"MPI [-]")
plt.ylabel(r"Efficiency [\%]")
plt.legend()
plt.tight_layout()
plt.savefig("MPI_Timing_Data_efficiency.pdf", bbox_inches='tight', dpi=300)
plt.close()