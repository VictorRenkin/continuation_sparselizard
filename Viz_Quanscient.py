import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# Définition globale des paramètres de police et de taille pour tous les graphiques
plt.rc('font', family='serif')  # Police avec empattements, comme Times
plt.rc('text', usetex=True)  # Utiliser LaTeX pour le texte dans les figures
plt.rcParams.update({
    'font.size': 18,         # Taille générale de la police
    'legend.fontsize': 28,   # Taille de la police des légendes
    'axes.labelsize': 25,    # Taille de la police des étiquettes des axes
    'axes.titlesize': 14,    # Taille de la police pour les titres des axes
    'xtick.labelsize': 14,   # Taille de la police pour les graduations de l'axe des X
    'ytick.labelsize': 14,   # Taille de la police pour les graduations de l'axe des Y
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

# df = pd.read_csv('data/Quascient_data/output_Pseudo-arcegth rapport mesh.csv')
# # Concaténation des vecteurs pour obtenir un seul vecteur frequency et forward
# forward_1, frequency_1 = extract_value_and_time(df, "Forward-1")
# forward_2, frequency_2 = extract_value_and_time(df, "Forward-2")
# forward_3, frequency_3 = extract_value_and_time(df, "Forward-3")

# # Application au vecteur fréquence
# frequency = smart_concatenate(frequency_1, frequency_2)
# frequency = smart_concatenate(frequency, frequency_3)

# # Et pour le vecteur forward correspondant
# forward = smart_concatenate(forward_1, forward_2)
# forward = smart_concatenate(forward, forward_3)
# bifurcation, frequency_bif = extract_value_and_time(df, "Bifurcation")


# # Tes données
# x_up = np.array([
#     852.5179856115108, 893.1654676258993, 925.8992805755396, 930.2158273381295,
#     934.5323741007194, 937.7697841726618, 941.726618705036, 949.6402877697842,
#     957.1942446043165, 966.1870503597122, 970.5035971223022, 974.4604316546763,
#     978.0575539568346, 981.6546762589928, 986.3309352517986, 989.568345323741,
#     993.5251798561151, 996.0431654676258, 998.5611510791367, 1001.4388489208634,
#     1003.9568345323742, 1006.8345323741007, 1008.9928057553957, 1010.7913669064748,
#     1012.9496402877699, 1015.1079136690647, 1017.2661870503598, 1020.5035971223023,
#     1023.3812949640287, 1026.978417266187, 1029.8561151079136, 1032.0143884892086,
#     1034.1726618705036, 1036.6906474820144, 1039.928057553957, 1042.4460431654677,
#     1045.68345323741, 1050, 1054.6762589928057, 1057.9136690647483,
#     1063.6690647482014, 1068.705035971223, 1072.6618705035971, 1077.3381294964029,
#     1083.093525179856, 1090.2877697841727, 1099.2805755395684, 1106.1151079136691,
#     1114.7482014388488, 1124.1007194244605, 1136.6906474820144
# ])

# y_up = np.array([
#     0.00019417475728155336, 0.0003155339805825242, 0.0004126213592233009,
#     0.0004368932038834951, 0.00046116504854368923, 0.0004854368932038834,
#     0.0005097087378640776, 0.0005582524271844659, 0.0006553398058252426,
#     0.0007766990291262134, 0.0008980582524271844, 0.000995145631067961,
#     0.0010922330097087377, 0.0011893203883495144, 0.0013592233009708736,
#     0.0015291262135922329, 0.0017718446601941746, 0.0020388349514563102,
#     0.0023543689320388345, 0.002694174757281553, 0.0030339805825242714,
#     0.003543689320388349, 0.003980582524271844, 0.004368932038834951,
#     0.004805825242718446, 0.005388349514563106, 0.005849514563106796,
#     0.006480582524271843, 0.007087378640776698, 0.007766990291262135,
#     0.008325242718446603, 0.00871359223300971, 0.009101941747572817,
#     0.009563106796116506, 0.010024271844660195, 0.01046116504854369,
#     0.010946601941747573, 0.01150485436893204, 0.012135922330097087,
#     0.012669902912621359, 0.01337378640776699, 0.014004854368932038,
#     0.014490291262135923, 0.015024271844660195, 0.015728155339805826,
#     0.016504854368932037, 0.01742718446601942, 0.018106796116504856,
#     0.018980582524271847, 0.01987864077669903, 0.020946601941747574
# ])

# # Division de x par 2*pi
# x_up = x_up / (2 * np.pi)



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
# plt.figure(figsize=(12, 8))  # Augmenter la taille de la figure pour une meilleure visibilité
# plt.scatter(frequency_bif, bifurcation, color=color_list[1], label=r'Bifurcation', marker='x', s=60, zorder=3)  # Augmenter la taille des marqueurs
# plt.plot(frequency, forward, color=color_list[0], linestyle='-', linewidth=2.5, label=r'Forward', zorder=1)  # Ligne continue plus épaisse
# plt.plot(x_up, y_up, color=color_list[2], linestyle='--', linewidth=1.5, label=r'Blahos 2022 (Up)', zorder=2)  # Ligne pointillée pour les données "up"
# plt.plot(x_down, y_down, color=color_list[2], linestyle='-.', linewidth=1.5, zorder=2)  # Ligne point-trait pour les données "down"
# plt.xlabel(r"Frequency [Hz]")
# plt.ylabel(r"Max displacement [m]")
# plt.xlim(np.min(frequency), np.max(frequency))  # Ajouter une marge pour mieux voir les extrémités
# plt.ylim(0, np.max(forward) * 1.2)  # Ajouter une marge pour mieux voir les données
# plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)  # Ajouter une grille pour une meilleure lisibilité
# plt.legend(loc='best', fontsize=14)  # Ajuster la taille de la légende
# plt.tight_layout()  # Optimiser l'agencement des éléments
# plt.savefig("NLFR_enhanced.pdf", bbox_inches='tight', dpi=300)

x2 = np.array([
    1.0171343873517782, 1.016679841897233, 1.0163043478260867, 1.0156916996047427,
    1.0151185770750986, 1.0144861660079048, 1.014071146245059, 1.0136758893280629,
    1.0133003952569166, 1.0129446640316202, 1.0125494071146242, 1.012154150197628,
    1.0118972332015808, 1.0116403162055334, 1.0115019762845847, 1.011363636363636,
    1.0112450592885374, 1.0111264822134385, 1.0110869565217389
]) * 24.8945/(np.pi * 2)
y2 = [
    0.17555178268251276, 0.17962648556876062, 0.18438030560271648, 0.19185059422750425,
    0.200679117147708, 0.2101867572156197, 0.2169779286926995, 0.22580645161290325,
    0.23531409168081496, 0.24278438030560273, 0.25500848896434636, 0.265874363327674,
    0.27741935483870966, 0.28828522920203736, 0.2977928692699491, 0.30662139219015283,
    0.3154499151103565, 0.3242784380305603, 0.33174872665534805
]

x = np.array([
    0.9916996047430832, 0.9922727272727274, 0.9930632411067195, 0.9939525691699606,
    0.9994861660079052, 1.0054743083003952, 1.0075296442687744, 1.0081422924901184,
    1.0084584980237152, 1.0087944664031618, 1.009169960474308, 1.0095256916996045,
    1.0099011857707507, 1.0102569169960471, 1.0106521739130432, 1.0110276679841894,
    1.0114031620553356, 1.0117193675889327, 1.0118379446640313
])* 24.8945/(np.pi * 2)
y = [
    0.24278438030560273, 0.2502546689303905, 0.26112054329371814, 0.2733446519524618,
    0.35755517826825123, 0.4519524617996604, 0.4825127334465195, 0.4899830220713073,
    0.4940577249575552, 0.498132427843803, 0.5022071307300509, 0.5076400679117148,
    0.5103565365025466, 0.5144312393887945, 0.5171477079796265, 0.5185059422750424,
    0.5191850594227504, 0.5157894736842106, 0.5117147707979626
]

df = pd.read_csv('data/Quascient_data/output_H0123456_cantivilateur.csv')
# Concaténation des vecteurs pour obtenir un seul vecteur frequency et forward
backward_1, frequency_1 = extract_value_and_time(df, "Backward-1")
backward_2, frequency_2 = extract_value_and_time(df, "Backward-2")
backward_3, frequency_3 = extract_value_and_time(df, "Backward-3")



# plt.figure(figsize=(12, 8))
# plt.plot(frequency_1, backward_1, color=color_list[0], label="Backward-1", marker='o', markersize=4, linewidth=1.5)
# plt.plot(frequency_2, backward_2, color=color_list[1], label="Backward-2", marker='o', markersize=4, linewidth=1.5)
# plt.plot(frequency_3, backward_3, color=color_list[2], label="Backward-3", marker='o', markersize=4, linewidth=1.5)
# plt.plot(x, y, color=color_list[3], markersize=4, linewidth=1.5)
# plt.plot(x2, y2, color=color_list[4], markersize=4, linewidth=1.5)
# plt.xlabel(r"Frequency [Hz]")
# plt.ylabel(r"Max displacement [m]")
# plt.legend(loc='best', fontsize=14)
# plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
# plt.tight_layout()
# plt.savefig("test_backward.pdf", bbox_inches='tight', dpi=300)
# # Application au vecteur fréquence

frequency = smart_concatenate(frequency_3, frequency_2)

# frequency = smart_concatenate(frequency, frequency_3)
# # Et pour le vecteur forward correspondant
backward = smart_concatenate(backward_3, backward_2)
backward = smart_concatenate(backward, backward_1)
frequency = smart_concatenate(frequency, frequency_1)
bifurcation, frequency_bif = extract_value_and_time(df, "Bifurcation")
plt.figure(figsize=(12, 8))  # Augmenter la taille de la figure pour une meilleure visibilité
plt.plot(frequency, backward, color=color_list[0], marker='o', markersize=4, label=r'Backward', linewidth=1.5, zorder=2)  # Augmenter la taille des marqueurs et l'épaisseur de la ligne
plt.scatter(frequency_bif, bifurcation, color=color_list[1], label=r'Bifurcation', marker='x', s=60, zorder=3)  # Augmenter la taille des marqueurs
plt.xlim(3.9,4.05)
plt.ylim(0.1, np.max(backward) * 1.1)  # Ajouter une marge pour mieux voir les données
plt.xlabel(r"Frequency [Hz]")
plt.ylabel(r"Max displacement [m]")
plt.legend(loc='best')  # Ajouter la légende
plt.tight_layout()  # Optimiser l'agencement des éléments
plt.savefig("NLFR_backward_cantivilateur.pdf", bbox_inches='tight', dpi=300)
# backward = smart_concatenate(backward, backward_3)
# bifurcation, frequency_bif = extract_value_and_time(df, "Bifurcation")
# plt.figure(figsize=(12, 8))  # Augmenter la taille de la figure pour une meilleure visibilité
# plt.plot(frequency, backward, color=color_list[0], marker='o', markersize=4, label=r'Backward', linewidth=1.5, zorder=2)  # Augmenter la taille des marqueurs et l'épaisseur de la ligne
# plt.xlabel(r"Frequency [Hz]")
# plt.ylabel(r"Max displacement [m]")
# plt.xlim(np.min(frequency), np.max(frequency))  # Ajouter une marge pour mieux voir les extrémités
# plt.ylim(0, np.max(backward) * 1.2)  # Ajouter une marge pour mieux voir les données
# plt.savefig("NLFR_backward_cantivilateur.pdf", bbox_inches='tight', dpi=300)