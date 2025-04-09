import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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


def real_time_plot_data_FRF(PATH, PATH_NEWTHON=None):
    """
    Plots the Frequency Response Function (FRF) using data from forward, backward, and prediction paths.

    Parameters
    ----------
    path : str
        Path where the resulting plot will be saved.
    csv_path_forward : str
        Path to the CSV file containing the forward path data (NLFR).
    csv_path_backward : str
        Path to the CSV file containing the backward path data. If empty, it is skipped.
    csv_pred_path : str, optional
        Path to the CSV file containing the prediction data. If not provided or empty, prediction is skipped.

    Returns
    -------
    None
        The function saves the plot at the specified path.
    """
    plt.figure(figsize=(8, 6))
    path_forward = PATH.get('PATH_STORE_DATA_FORWARD')
    if path_forward:
        try:
            df_forward = pd.read_csv(path_forward)
            if not df_forward.empty:
                plt.plot(df_forward['freq'], df_forward['u'], color=color_list[0], marker='o', label=r'NLFR forward')
                last_u = df_forward['u'].iloc[-1]
                last_f = df_forward['freq'].iloc[-1]
        except FileNotFoundError:
            pass
    else:
        pass
    # Bloc backward
    path_backward = PATH.get('PATH_STORE_DATA_DOWNWARD')
    if path_backward:
        try:
            df_backward = pd.read_csv(path_backward)
            if not df_backward.empty:
                plt.plot(df_backward['freq'], df_backward['u'], color=color_list[1], marker='o', label=r'NLFR backward')
                last_u = df_backward['u'].iloc[-1]
                last_f = df_backward['freq'].iloc[-1]
        except FileNotFoundError:
            pass
    else:
        pass

    # Bloc predictor
    path_predictor = PATH.get('PATH_STORE_PREDICTOR')
    if path_predictor:
        try:
            df_pred = pd.read_csv(path_predictor)
            if not df_pred.empty:
                u_pred = df_pred['u'].iloc[-1]
                f_pred = df_pred['freq'].iloc[-1]
                plt.plot([last_f, f_pred], [last_u, u_pred], color=color_list[2], linestyle='--', marker='o', label=r'Prediction step')
        except FileNotFoundError:
            pass
    else:
        pass

        
    if PATH_NEWTHON:
        try:
            df_newthon = pd.read_csv(PATH_NEWTHON)
            plt.scatter(df_newthon['freq'], df_newthon['u'], facecolors='none', edgecolors=color_list[3], label = r'Newton iteration')
        except FileNotFoundError:
            print(f"Prediction file not found: {PATH['PATH_STORE_DATA_FORWARD']}")


    plt.xlabel(r"Frequency [Hz]")
    plt.ylabel(r"Amplitude [m]")
    plt.legend()
    plt.tight_layout()
    plt.savefig(PATH["PATH_FIGURE"], format='pdf', bbox_inches='tight', transparent=True)
    plt.close()

def viz_forward_and_backward(PATH_FORWARD, PATH_DOWNWARD, PATH_FIGURE) : 
    df_forward = pd.read_csv(PATH_FORWARD)
    df_downward = pd.read_csv(PATH_DOWNWARD)

    plt.plot(df_forward['freq'], df_forward['u'], color=color_list[0], marker='o', label=r'NLFR forward')  
    plt.plot(df_downward['freq'], df_downward['u'], color=color_list[1], marker='o', label=r'NLFR downward')  

    plt.xlabel(r"Frequency [Hz]")
    plt.ylabel(r"Amplitude [m]")
    plt.legend()
    plt.tight_layout()
    plt.savefig(PATH_FIGURE, format='pdf', bbox_inches='tight', transparent=True)
    plt.close()




