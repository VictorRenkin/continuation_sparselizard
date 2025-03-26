import os
import pandas as pd


def create_doc_csv(path):
    """
    Creates a new CSV file at the specified path. If the file already exists, it is deleted.
    If the directory does not exist, it is created before generating the CSV file.
    
    Parameters:
    ----------
    path (str): 
        The path to the CSV file.
    """
    directory = os.path.dirname(path)
    
    if directory and not os.path.exists(directory):
        os.makedirs(directory)  # Create the directory if it doesn't exist

    if os.path.exists(path):
        os.remove(path)  # Delete existing file

    df = pd.DataFrame(columns=["u", "freq"])
    df.to_csv(path, index=False)  # Create the new CSV f

def add_data_to_csv(u_max, freq, path):
    """
    Appends a new row with 'u_max' and 'freq' values to the specified CSV file.

    Parameters:
    u_max (float): 
        The maximum value of 'u' to be recorded.
    freq (float): 
        The frequency value to be recorded.
    path (str): 
        The path to the CSV file.
    """
    df = pd.DataFrame({"freq": [freq], "u": [u_max]})
    df.to_csv(path, mode='a', header=False, index=False)

def remove_last_row_from_csv(path):
    """
    Removes the last row from the CSV file at the given path.
    
    Parameters:
    path (str): 
        The path to the CSV file.
    """
    df = pd.read_csv(path)
    df = df.iloc[:-1]  # enlève la dernière ligne
    df.to_csv(path, index=False)