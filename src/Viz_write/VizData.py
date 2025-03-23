from import_extension.imports import *

# Définition globale des paramètres de police et de taille pour tous les graphiques
plt.rc('font', family='serif')  # Police avec empattements, comme Times
plt.rc('text', usetex=True)  # Utiliser LaTeX pour le texte dans les figures
plt.rcParams.update({
    'font.size': 18,         # Taille générale de la police
    'legend.fontsize': 18,   # Taille de la police des légendes
    'axes.labelsize': 18,    # Taille de la police des étiquettes des axes
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

def viz_convergence(path="../figures/"):
    # Données à visualiser

    corrected_frequency_data = [
    [183.759, 183.759, 187.276, 187.276],
    [162.764, 162.764, 163.923, 163.923],
    [160.446, 160.446, 160.767, 160.767],
    [160.122, 160.122, 160.333, 160.333],
    [159.587, 159.587, 159.629, 159.629],
    [159.478, 159.478, 159.481, 159.481]]

    number_node = [75, 165, 315, 405, 1265, 21825]
    corrected_relative_sums = []

    for i in range(len(corrected_frequency_data)-1):
        base_frequency = corrected_frequency_data[-1]  # Référence
        relative_diff_sum = np.abs(corrected_frequency_data[i][0] - base_frequency[0]) /  base_frequency[0]  * 100
        corrected_relative_sums.append(relative_diff_sum)

    plt.figure(figsize=(8, 5))
    plt.xlim(number_node[0], number_node[-2])
    plt.plot(number_node[:-1], corrected_relative_sums, marker='o', linestyle='-', color=color_list[0])
    plt.hlines(3,number_node[0],number_node[-2],linestyles='-.',color=color_list[1])
    plt.text(number_node[-2]-750, 3.3,r'Acceptable range',color='black')
    plt.xscale("log")
    plt.xlabel(r"Number of node [-]")
    plt.ylabel(r"Relative error [\%]")
    plt.savefig(path+"convergence.pdf", bbox_inches='tight', dpi=300)
    plt.close()

def viz_NLFR(freq, amplitude, path = '../figures/') :
    plt.figure()
    plt.plot(freq, amplitude, color = color_list[0], marker='o')
    plt.xlabel(r"Frequency [Hz]")
    plt.ylabel(r"Amplitude [m]")
    plt.savefig(path+"NLFR.pdf", bbox_inches='tight', dpi=300)
    plt.close()

def real_time_plot_data_FRF(path, csv_path):
    df = pd.read_csv(csv_path)
    plt.figure(figsize=(8, 6))
    plt.plot(df["u"], df["freq"], color = color_list[0], marker='o')
    plt.xlabel(r"Frequency [Hz]")
    plt.ylabel(r"Amplitude [m]")
    plt.savefig(path)  # Save the updated figure
    plt.close()  # Close the figure to free memory
    
