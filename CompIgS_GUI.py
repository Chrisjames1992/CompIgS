#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QLineEdit, QFileDialog, QProgressBar, QWidget, QFrame, QMessageBox
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt5.QtGui import QFont, QIcon, QPalette, QColor, QLinearGradient, QBrush
import seaborn as sns
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import peptides as py
from peptides import Peptide
from collections import Counter
from io import StringIO
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtWidgets import QDialog




class AnalysisThread(QThread):
    progress_updated = pyqtSignal(int)
    analysis_completed = pyqtSignal(str)
    error_occurred = pyqtSignal(str)

    def __init__(self, file_paths):
        super().__init__()
        self.file_paths = file_paths

    def run(self):
        try:
            import pandas as pd
            import os
            from Bio import SeqIO
        # Simulate file processing with progress updates
            for i in range(1, 101):
                self.msleep(20)  # Simulate processing delay
                self.progress_updated.emit(i)
            

            # Example of loading files into DataFrames (replace with actual analysis logic)
            stats_IgE = pd.read_csv(self.file_paths[0], sep='\t', low_memory=False)
            stats_IgG1 = pd.read_csv(self.file_paths[1], sep='\t', low_memory=False)
            stats_IgE_mut = pd.read_csv(self.file_paths[2], sep='\t', low_memory=False)
            stats_IgG1_mut = pd.read_csv(self.file_paths[3], sep='\t', low_memory=False)
            stats_IgE_Aminotab = pd.read_csv(self.file_paths[4], sep='\t', low_memory=False)
            stats_IgG1_Aminotab = pd.read_csv(self.file_paths[5], sep='\t', low_memory=False)
            stats_IgE_Aminotab_change = pd.read_csv(self.file_paths[6], sep='\t', low_memory=False)
            stats_IgG1_Aminotab_change = pd.read_csv(self.file_paths[7], sep='\t', low_memory=False)
            fasta_file_path = self.file_paths[8]
            sequences = list(SeqIO.parse(fasta_file_path, "fasta"))
            
            # Analysis logic starts here
            stats_IgE_productive = stats_IgE[(stats_IgE['functionality'] == "productive") &
                                             (stats_IgE['onecopy'] >= 2) &
                                             ~(stats_IgE['vgene'] == '') & ~(stats_IgE['vgene'].isna()) &
                                             ~(stats_IgE['dgene'] == '') & ~(stats_IgE['dgene'].isna()) &
                                             ~(stats_IgE['jgene'] == '') & ~(stats_IgE['jgene'].isna())]
    
             # Create a copy of the productive IgE DataFrame
            stats_IgE_productive_VDJ = stats_IgE_productive.copy()
            
            # Check if the input DataFrame is empty before proceeding
            if not stats_IgE_productive.empty:
                # Concatenate vgene, dgene, and jgene to create the 'IgE_VDJ' column
                stats_IgE_productive_VDJ = stats_IgE_productive.copy()
                stats_IgE_productive_VDJ['IgE_VDJ'] = stats_IgE_productive['vgene'].astype(str) + \
                                                      stats_IgE_productive['dgene'].astype(str) + \
                                                      stats_IgE_productive['jgene'].astype(str)
            
                # Remove spaces in the 'IgE_VDJ' column
                stats_IgE_productive_VDJ['IgE_VDJ'] = stats_IgE_productive_VDJ['IgE_VDJ'].str.replace(' ', '', regex=False)
            
                # Calculate the number of unique IgE clones
                unique_strings_IgE_productive_VDJ = stats_IgE_productive_VDJ['IgE_VDJ'].nunique(dropna=True)
            
                # Check if unique clones were found
                if unique_strings_IgE_productive_VDJ > 0:
                    print(f'Number of IgE clones: {unique_strings_IgE_productive_VDJ}')
                else:
                    print("No IgE_VDJ clones found.")
            else:
                print("No entries found in the input DataFrame 'stats_IgE_productive'. Continuing analysis with zero clones.")
                unique_strings_IgE_productive_VDJ = 0


    
            # Filter productive IgG1 sequences with additional criteria
            stats_IgG1_productive = stats_IgG1[(stats_IgG1['functionality'] == "productive") &
                                               (stats_IgG1['onecopy'] >= 2) &
                                               ~(stats_IgG1['vgene'] == '') & ~(stats_IgG1['vgene'].isna()) &
                                               ~(stats_IgG1['dgene'] == '') & ~(stats_IgG1['dgene'].isna()) &
                                               ~(stats_IgG1['jgene'] == '') & ~(stats_IgG1['jgene'].isna())]
            
            # Check if the filtered DataFrame is empty
            if not stats_IgG1_productive.empty:
                # Create a copy and generate the 'IgG1_VDJ' column by concatenating vgene, dgene, and jgene
                stats_IgG1_productive_VDJ = stats_IgG1_productive.copy()
                stats_IgG1_productive_VDJ['IgG1_VDJ'] = stats_IgG1_productive['vgene'].astype(str) + \
                                                        stats_IgG1_productive['dgene'].astype(str) + \
                                                        stats_IgG1_productive['jgene'].astype(str)
            
                # Remove spaces from the 'IgG1_VDJ' column
                stats_IgG1_productive_VDJ['IgG1_VDJ'] = stats_IgG1_productive_VDJ['IgG1_VDJ'].str.replace(' ', '', regex=False)
            
                # Calculate the number of unique IgG1 clones
                unique_strings_IgG1_productive_VDJ = stats_IgG1_productive_VDJ['IgG1_VDJ'].nunique(dropna=True)
                print(f'Number of IgG1 clones: {unique_strings_IgG1_productive_VDJ}')
            else:
                # Handle the case where no productive sequences are found
                print("No entries found in the input DataFrame 'stats_IgG1_productive'. Continuing analysis with zero clones.")
                unique_strings_IgG1_productive_VDJ = 0

    
            # Emit signal after successful analysis
            self.analysis_completed.emit("Analysis completed successfully!")
    
            # Make a copy of the stats_IgE_productive_VDJ
            stats_IgE_productive_VDJ_compare = stats_IgE_productive_VDJ.copy()
            stats_IgG1_productive_VDJ_compare = stats_IgG1_productive_VDJ.copy()
            
            # Compare the 'IgE VDJ and IgG1 VDJ' values between the two DataFrames and add the result as a new column in each dataframe
            stats_IgE_productive_VDJ_compare['Match_in_IgG1_VDJ'] = stats_IgE_productive_VDJ_compare['IgE_VDJ'].apply(lambda x: x in set(stats_IgG1_productive_VDJ_compare['IgG1_VDJ']))
            stats_IgG1_productive_VDJ_compare['Match_in_IgE_VDJ'] = stats_IgG1_productive_VDJ_compare['IgG1_VDJ'].apply(lambda x: x in set(stats_IgE_productive_VDJ_compare['IgE_VDJ']))
            
            #Filter out shared IgE clones,number of clones and percentage of shared IgE clones
            shared_IgE = stats_IgE_productive_VDJ_compare[stats_IgE_productive_VDJ_compare['Match_in_IgG1_VDJ']]
            unique_strings_shared_IgE = shared_IgE['IgE_VDJ'].nunique(dropna=True)
            print(f'Number of shared IgE clones: {unique_strings_shared_IgE}')
            percentage_shared_IgE_clones = (unique_strings_shared_IgE / unique_strings_IgE_productive_VDJ) * 100
            print(f'percentage of shared IgE clones: {percentage_shared_IgE_clones}')
            
            #Filter out Unique IgE clones,number of clones and percentage of Unique IgE clones
            Unique_IgE = stats_IgE_productive_VDJ_compare[~stats_IgE_productive_VDJ_compare['Match_in_IgG1_VDJ']]
            unique_strings_unique_IgE = Unique_IgE['IgE_VDJ'].nunique(dropna=True)
            print(f'Number of unique IgE clones: {unique_strings_unique_IgE}')
            percentage_unique_IgE_clones = (unique_strings_unique_IgE / unique_strings_IgE_productive_VDJ) * 100
            print(f'percentage of unique IgE clones: {percentage_unique_IgE_clones}')
            
            #Filter out shared IgG1 clones,number of clones and percentage of shared IgG1 clones
            shared_IgG1 = stats_IgG1_productive_VDJ_compare[stats_IgG1_productive_VDJ_compare['Match_in_IgE_VDJ']]
            unique_strings_shared_IgG1 = shared_IgG1['IgG1_VDJ'].nunique(dropna=True)
            print(f'Number of shared IgG1 clones: {unique_strings_shared_IgG1}')
            percentage_shared_IgG1_clones = (unique_strings_shared_IgG1 / unique_strings_IgG1_productive_VDJ) * 100
            print(f'percentage of shared IgG1 clones: {percentage_shared_IgG1_clones}')
            
            #Filter out Unique IgG1 clones,number of clones and percentage of unique IgG1 clones
            Unique_IgG1 = stats_IgG1_productive_VDJ_compare[~stats_IgG1_productive_VDJ_compare['Match_in_IgE_VDJ']]
            unique_strings_unique_IgG1 = Unique_IgG1['IgG1_VDJ'].nunique(dropna=True)
            print(f'Number of unique IgG1 clones: {unique_strings_unique_IgG1}')
            percentage_unique_IgG1_clones = (unique_strings_unique_IgG1 / unique_strings_IgG1_productive_VDJ) * 100
            print(f'percentage of unique IgG1 clones: {percentage_unique_IgG1_clones}')

            
                
            #Plot a graph of number of shared and unique IgE and IgG1 clones
            # Define the data
            types = ['IgE','IgG1', 'IgE','IgE','IgG1','IgG1']
            categories = ['Total IgE', 'Total IgG1', 'Shared IgE', 'Unique IgE', 'Shared IgG1', 'Unique IgG1']
            counts = [unique_strings_IgE_productive_VDJ, unique_strings_IgG1_productive_VDJ, unique_strings_shared_IgE, unique_strings_unique_IgE,unique_strings_shared_IgG1, unique_strings_unique_IgG1]
            
            # Create a DataFrame
            Graph1 = pd.DataFrame({
            'Type': types,
            'Category': categories,
            'Count': counts
            })
            # Set Seaborn style for scientific plots
            sns.set(style="whitegrid")
            
            # Filter data for IgE and IgG1
            df_IgE = Graph1[Graph1['Type'] == 'IgE']
            df_IgG1 = Graph1[Graph1['Type'] == 'IgG1']
            
            # Calculate dynamic y-axis limits
            max_IgE = df_IgE['Count'].max()
            max_IgG1 = df_IgG1['Count'].max()
            
            ylim_IgE = max_IgE + 10  # Slightly above the max value
            ylim_IgG1 = max_IgG1 + 100  # Slightly above the max value
            
            # Generate dynamic color palette based on the number of unique categories
            unique_categories_IgE = df_IgE['Category'].nunique()
            unique_categories_IgG1 = df_IgG1['Category'].nunique()
            
            # Get n unique colors for each plot
            palette_IgE = sns.color_palette("Set2", unique_categories_IgE)
            palette_IgG1 = sns.color_palette("Set2", unique_categories_IgG1)
            
            # Create two subplots
            fig1, axes = plt.subplots(1, 2, figsize=(12, 6))
            
            # Plot for IgE with dynamically generated color palette
            sns.barplot(
            data=df_IgE,
            x='Category',
            y='Count',
            hue='Category',
            palette=palette_IgE,
            dodge=False,
            ax=axes[0]
            )
            axes[0].set_title('Graph A: IgE', fontsize=14)
            axes[0].set_xlabel('', fontsize=12)
            axes[0].set_ylabel('Number of Clones', fontsize=12)
            axes[0].set_ylim(0, ylim_IgE)  # Set dynamic y-axis limit for IgE
            axes[0].legend([], [], frameon=False)  # Remove legend
            for container in axes[0].containers:
                axes[0].bar_label(container, fmt='%d', label_type='edge', fontsize=10)
            
            # Plot for IgG1 with dynamically generated color palette
            sns.barplot(
            data=df_IgG1,
            x='Category',
            y='Count',
            hue='Category',
            palette=palette_IgG1,
            dodge=False,
            ax=axes[1]
            )
            axes[1].set_title('Graph B: IgG1', fontsize=14)
            axes[1].set_xlabel('', fontsize=12)
            axes[1].set_ylabel('')  # No ylabel for the second plot
            axes[1].set_ylim(0, ylim_IgG1)  # Set dynamic y-axis limit for IgG1
            axes[1].legend([], [], frameon=False)  # Remove legend
            for container in axes[1].containers:
                axes[1].bar_label(container, fmt='%d', label_type='edge', fontsize=10)
    
            
            # Adjust layout
            plt.tight_layout()  
            
            # Determine the Downloads folder of the current user
            downloads_folder = os.path.join(os.path.expanduser("~"), "Downloads")
            
            # Create a folder for the analysis plots
            plots_folder = os.path.join(downloads_folder, "BCR_Analysis_Plots","Shared IgE")
            os.makedirs(plots_folder, exist_ok=True)  # Create all intermediate directories if they don't exist
            
            # Define the plot path
            plot_path = os.path.join(plots_folder, "Number of shared and unique clones.png")
            
            # Save the plot
            plt.savefig(plot_path, dpi=1200, bbox_inches='tight')
            
            
            
            # Sort out shared IgE for VDJ rearrangement with two V gene call
            # Filter rows where 'IgE_VDJ' contains 'or', ignoring NaNs
            shared_IgE_two_vgene = shared_IgE.loc[
                shared_IgE['IgE_VDJ'].str.contains('or', case=False, na=False)
            ].copy()
            
            # Apply the transformation to the specific column 'IgE_VDJ'
            shared_IgE_two_vgene['IgE_VDJ'] = shared_IgE_two_vgene['IgE_VDJ'].apply(
                lambda x: x.split('or', 1)[-1].strip() if isinstance(x, str) and 'or' in x else x
            )
            
            # Find IgE VDJ from shared_IgE_two_vgene that are still represented in shared IgE
            shared_IgE_two_vgene_1 = shared_IgE_two_vgene.copy()
            shared_IgE_two_vgene_1['Match_in_Shared_IgE_VDJ'] = shared_IgE_two_vgene_1['IgE_VDJ'].apply(lambda x: x in set(shared_IgE['IgE_VDJ']))
            
            # Check if there are matches before proceeding
            if shared_IgE_two_vgene_1['Match_in_Shared_IgE_VDJ'].any():
                # Filter out the processed two v gene IgE VDJ that are in shared IgE
                shared_IgE_two_vgene_2 = shared_IgE_two_vgene_1[shared_IgE_two_vgene_1['Match_in_Shared_IgE_VDJ']]
                
                # Drop 'Match_in_Shared_IgE_VDJ' column
                shared_IgE_two_vgene_3 = shared_IgE_two_vgene_2.drop(columns=['Match_in_Shared_IgE_VDJ'])
                
                # Merge the filtered out two v gene IgE VDJ with shared IgE to create a new dataframe
                shared_IgE_main = pd.concat([shared_IgE, shared_IgE_two_vgene_3], ignore_index=True)
            else:
                print("No matches found in shared IgE VDJ.")
                shared_IgE_main = shared_IgE.copy()  # Continue with the original shared_IgE if no matches

            
            
            # Sort out shared IgG1 for VDJ rearrangement with two V gene call
            # Filter rows where 'IgG1_VDJ' contains 'or', ignoring NaNs
            shared_IgG1_two_vgene = shared_IgG1.loc[
                shared_IgG1['IgG1_VDJ'].str.contains('or', case=False, na=False)
            ].copy()
            
            # Apply the transformation to the specific column 'IgG1_VDJ'
            shared_IgG1_two_vgene['IgG1_VDJ'] = shared_IgG1_two_vgene['IgG1_VDJ'].apply(
                lambda x: x.split('or', 1)[-1].strip() if isinstance(x, str) and 'or' in x else x
            )
            
            # Find IgG1 VDJ from shared_IgG1_two_vgene that are still represented in shared IgG1
            shared_IgG1_two_vgene_1 = shared_IgG1_two_vgene.copy()
            shared_IgG1_two_vgene_1['Match_in_Shared_IgG1_VDJ'] = shared_IgG1_two_vgene_1['IgG1_VDJ'].apply(lambda x: x in set(shared_IgG1['IgG1_VDJ']))
            
            # Check if there are matches before proceeding
            if shared_IgG1_two_vgene_1['Match_in_Shared_IgG1_VDJ'].any():
                # Filter out the processed two v gene IgG1 VDJ that are in shared IgG1
                shared_IgG1_two_vgene_2 = shared_IgG1_two_vgene_1[shared_IgG1_two_vgene_1['Match_in_Shared_IgG1_VDJ']]
                
                # Drop 'Match_in_Shared_IgG1_VDJ' column
                shared_IgG1_two_vgene_3 = shared_IgG1_two_vgene_2.drop(columns=['Match_in_Shared_IgG1_VDJ'])
                
                # Merge the filtered out two v gene IgG1 VDJ with shared IgG1 to create a new dataframe
                shared_IgG1_main = pd.concat([shared_IgG1, shared_IgG1_two_vgene_3], ignore_index=True)
            else:
                print("No matches found in shared IgG1 VDJ.")
                shared_IgG1_main = shared_IgG1.copy()  # Continue with the original shared_IgG1 if no matches

            
            #Extract the cd3raa, onecopy, IgE_VDJ columns
            shared_IgE_main_1 = shared_IgE_main[['cdr3aa', 'onecopy', 'IgE_VDJ']]
            
            # Group by 'IgE_VDJ' and sum the 'onecopy' values
            subtotal_shared_IgE = shared_IgE_main_1.groupby('IgE_VDJ', as_index=False)['onecopy'].sum()
            
            # Rename the column for clarity
            subtotal_shared_IgE.rename(columns={'onecopy': 'subtotal_onecopy'}, inplace=True)
            
            #Extract the cd3raa, onecopy, IgG1_VDJ columns
            shared_IgG1_main_1 = shared_IgG1_main[['cdr3aa', 'onecopy', 'IgG1_VDJ']]
            
            # Group by 'IgG1_VDJ' and sum the 'onecopy' values
            subtotal_shared_IgG1 = shared_IgG1_main_1.groupby('IgG1_VDJ', as_index=False)['onecopy'].sum()
            
            # Rename the column for clarity
            subtotal_shared_IgG1.rename(columns={'onecopy': 'subtotal_onecopy'}, inplace=True)
            
            # Calculate the ratio of logIgE/IgG1
            # Ensure both DataFrames have the same number of rows or handle them appropriately
            if len(subtotal_shared_IgE) == len(subtotal_shared_IgG1):
                ratio1 = subtotal_shared_IgE['subtotal_onecopy'] / subtotal_shared_IgG1['subtotal_onecopy']
                ratio2 = subtotal_shared_IgG1['subtotal_onecopy'] / subtotal_shared_IgE['subtotal_onecopy']
                # Add the ratio as a new column in both DataFrames
                subtotal_shared_IgE['IgE_IgG1_ratio'] = ratio1
                subtotal_shared_IgG1['IgG1_IgE_ratio'] = ratio2
                print("DataFrames must have the same number of rows to calculate the ratio.")
            
            # Log transform the IgE and IgG1 subtotal_onecopy
            
            subtotal_shared_IgE_log = subtotal_shared_IgE.copy()
            subtotal_shared_IgG1_log = subtotal_shared_IgG1.copy()
            subtotal_shared_IgE_log['log10_IgE_IgG1_onecopy'] = np.log10(subtotal_shared_IgE_log['IgE_IgG1_ratio'])
            subtotal_shared_IgG1_log['log10_IgG1_IgE_onecopy'] = np.log10(subtotal_shared_IgG1_log['IgG1_IgE_ratio'])
            
            # Sort the DataFrame in descending order based on log10_IgE_IgG1_onecopy
            subtotal_shared_IgE_log_sorted = subtotal_shared_IgE_log.sort_values(by='log10_IgE_IgG1_onecopy', ascending=False)
        
            # Set Seaborn style
            sns.set(style="whitegrid")


            
                
            # Create a color map for positive and negative values
            colors = subtotal_shared_IgE_log_sorted['log10_IgE_IgG1_onecopy'].apply(lambda y: 'green' if y >= 0 else 'red')
            
            # Plot the bar graph
            fig2 = plt.figure(figsize=(12, 6))
            bars = plt.bar(subtotal_shared_IgE_log_sorted['IgE_VDJ'],
            subtotal_shared_IgE_log_sorted['log10_IgE_IgG1_onecopy'],
            color=colors)
            
            # Add labels
            plt.xlabel('VH Gene', fontsize=12)
            plt.ylabel('Log10(IgE/IgG1)', fontsize=12)
            plt.title('Ratio of IgE/IgG1 Copy Number', fontsize=14)
            
            # Rotate x-axis labels for better readability
            plt.xticks(rotation=90, ha='center', fontsize=4)
            
            # Add a legend
            legend_labels = ['IgE-biased clones', 'IgG1-biased clones']
            handles = [plt.Rectangle((0, 0), 1, 1, color='green'), plt.Rectangle((0, 0), 1, 1, color='red')]
            plt.legend(handles, legend_labels, title="", fontsize=10)
            
        
            # Adjust layout
            plt.tight_layout() 
            
            downloads_folder = os.path.join(os.path.expanduser("~"), "Downloads")
            
            plots_folder = os.path.join(downloads_folder, "BCR_Analysis_Plots")
            if not os.path.exists(plots_folder):
                os.makedirs(plots_folder)  # Create the folder if it doesn't exist
                
            plot_path = os.path.join(plots_folder, "Shared IgE/IgG1 copy number ratio.png")
            plt.savefig(plot_path, dpi=1200, bbox_inches='tight')
            
            # Show the plot
                
                
            # Plot graph for cloal distribution of shared IgE and shared IgG1 copy number
            # Sort the counts in ascending order
            subtotal_shared_IgE_log_copy_sort = shared_IgE_main_1.sort_values(by='onecopy', ascending=True).reset_index(drop=True)
            subtotal_shared_IgG1_log_copy_sort = shared_IgG1_main_1.sort_values(by='onecopy', ascending=True).reset_index(drop=True)
            
            # Calculate the mean and median values for both dataframes
            mean1 = subtotal_shared_IgE_log_copy_sort['onecopy'].mean()
            mean2 = subtotal_shared_IgG1_log_copy_sort['onecopy'].mean()
            
            
            # Set Seaborn style
            sns.set(style="whitegrid")
            
            # Create the plot with two rows and two columns of subplots (original plots and zoomed-in plots)
            fig3, axes = plt.subplots(2, 2, figsize=(14, 12))
            
            # Plot for DataFrame 1 (Shared IgE) - Original
            sns.scatterplot(data=subtotal_shared_IgE_log_copy_sort, x=subtotal_shared_IgE_log_copy_sort.index, y='onecopy', color='blue', s=100, ax=axes[0, 0], label='Shared IgE', edgecolor='w')
            axes[0, 0].axhline(mean1, color='orange', linestyle='--', label=f'IgE Mean: {mean1:.2f}')
            axes[0, 0].axhline(mean2, color='green', linestyle='--', label=f'IgG1 Mean: {mean2:.2f}')
            axes[0, 0].set_title('Clonal Size Distribution - Shared IgE', fontsize=14)
            axes[0, 0].set_xlabel('Clonal Index', fontsize=12)
            axes[0, 0].set_ylabel('Clonal Size', fontsize=12)
            axes[0, 0].legend()
            
            # Plot for DataFrame 2 (Shared IgG1) - Original
            sns.scatterplot(data=subtotal_shared_IgG1_log_copy_sort, x=subtotal_shared_IgG1_log_copy_sort.index, y='onecopy', color='red', s=100, ax=axes[0, 1], label='Shared IgG1', edgecolor='w')
            axes[0, 1].axhline(mean2, color='green', linestyle='--', label=f'IgG1 Mean: {mean2:.2f}')
            axes[0, 1].axhline(mean1, color='orange', linestyle='--', label=f'IgE Mean: {mean1:.2f}')
            axes[0, 1].set_title('Clonal Size Distribution - Shared IgG1', fontsize=14)
            axes[0, 1].set_xlabel('Clonal Index', fontsize=12)
            axes[0, 1].set_ylabel('Clonal Size', fontsize=12)
            axes[0, 1].legend()
            
            # Define a range around the mean for zooming in (e.g., +-10 around the mean)
            zoom_range_1 = (mean2 - 20, mean2 + 20)
            zoom_range_2 = (mean2 - 20, mean2 + 20)
            
            # Filter DataFrames to show values around the mean
            shared_IgE_zoom = subtotal_shared_IgE_log_copy_sort[(subtotal_shared_IgE_log_copy_sort['onecopy'] >= zoom_range_1[0]) & (subtotal_shared_IgE_log_copy_sort['onecopy'] <= zoom_range_1[1])]
            shared_IgG1_zoom = subtotal_shared_IgG1_log_copy_sort[(subtotal_shared_IgG1_log_copy_sort['onecopy'] >= zoom_range_2[0]) & (subtotal_shared_IgG1_log_copy_sort['onecopy'] <= zoom_range_2[1])]
            
            # Plot for DataFrame 1 (Shared IgE) - Zoomed-in around the mean
            sns.scatterplot(data=shared_IgE_zoom, x=shared_IgE_zoom.index, y='onecopy', color='blue', s=100, ax=axes[1, 0], label='Shared IgE', edgecolor='w')
            axes[1, 0].axhline(mean2, color='green', linestyle='--', label=f'IgG1 Mean: {mean2:.2f}')
            axes[1, 0].axhline(mean1, color='orange', linestyle='--', label=f'IgE Mean: {mean1:.2f}')
            axes[1, 0].set_title('Zoomed Clonal Size Distribution - Shared IgE', fontsize=14)
            axes[1, 0].set_xlabel('Clonal Index', fontsize=12)
            axes[1, 0].set_ylabel('Clonal Size', fontsize=12)
            axes[1, 0].legend()
            
            # Plot for DataFrame 2 (Shared IgG1) - Zoomed-in around the mean
            sns.scatterplot(data=shared_IgG1_zoom, x=shared_IgG1_zoom.index, y='onecopy', color='red', s=100, ax=axes[1, 1], label='Shared IgG1', edgecolor='w')
            axes[1, 1].axhline(mean2, color='green', linestyle='--', label=f'IgG1 Mean: {mean2:.2f}')
            axes[1, 1].axhline(mean1, color='orange', linestyle='--', label=f'IgE Mean: {mean1:.2f}')
            axes[1, 1].set_title('Zoomed Clonal Size Distribution - Shared IgG1', fontsize=14)
            axes[1, 1].set_xlabel('Clonal Index', fontsize=12)
            axes[1, 1].set_ylabel('Clonal Size', fontsize=12)
            axes[1, 1].legend()
            
            # Adjust layout to avoid overlap
            plt.tight_layout()
            
            downloads_folder = os.path.join(os.path.expanduser("~"), "Downloads")
            
            plots_folder = os.path.join(downloads_folder, "BCR_Analysis_Plots")
            if not os.path.exists(plots_folder):
                os.makedirs(plots_folder)  # Create the folder if it doesn't exist
                
            plot_path = os.path.join(plots_folder, "Clonal size distribution of Shared IgE and IgG1.png")
            plt.savefig(plot_path, dpi=1200, bbox_inches='tight')
            
            # Show the plot
            
            
            
            
            #select IgE-biased and IgG1-biased clones
            
            #select IgE-biased clones
            # Create a copy of the dataframe
            IgE_biased_clones = subtotal_shared_IgE_log.copy()
            
            # Filter rows where 'log10_IgE_IgG1_onecopy' is greater than or equal to 0.301
            IgE_biased_clones = subtotal_shared_IgE_log[subtotal_shared_IgE_log['log10_IgE_IgG1_onecopy'] >= 0.301]
            
            #Compare the VDJ of the IgE-biased clones to the shared IgE VDJ
            IgE_IgG1_biased_clones_VDJ_compare = shared_IgE_main.copy()
            IgE_IgG1_biased_clones_VDJ_compare['Match_in_IgE_biased_VDJ'] = IgE_IgG1_biased_clones_VDJ_compare['IgE_VDJ'].apply(lambda x: x in set(IgE_biased_clones['IgE_VDJ']))
            
            #Sort out the IgE_biased clones details cd3raa,seqid......
            IgE_biased_clones_VDJ_compare = IgE_IgG1_biased_clones_VDJ_compare [IgE_IgG1_biased_clones_VDJ_compare['Match_in_IgE_biased_VDJ']]
            
            #select IgG1-biased clones
            IgG1_biased_clones = subtotal_shared_IgG1_log.copy()
            
            # Filter rows where 'log10_IgE_IgG1_onecopy' is greater than or equal to 0.301
            IgG1_biased_clones = subtotal_shared_IgG1_log[subtotal_shared_IgG1_log['log10_IgG1_IgE_onecopy'] >= 0.301]
            
            #Compare the VDJ of the IgG1-biased clones to the shared IgE VDJ
            IgE_IgG1_biased_clones_VDJ_compare['Match_in_IgG1_biased_VDJ'] = IgE_IgG1_biased_clones_VDJ_compare['IgE_VDJ'].apply(lambda x: x in set(IgG1_biased_clones['IgG1_VDJ']))
            
            #Sort out the IgG1_biased clones details cd3raa,seqid........
            IgG1_biased_clones_VDJ_compare = IgE_IgG1_biased_clones_VDJ_compare [IgE_IgG1_biased_clones_VDJ_compare['Match_in_IgG1_biased_VDJ']]
            
            # Cleanup IgE mutation table, remove Non-productive clones and clean up the number of mutations,CDR3 non-silent mutations and V region number of nucleotides
            # Filter for productive rows
            stats_IgE_mut_1 = stats_IgE_mut[stats_IgE_mut['V-DOMAIN Functionality'] == "productive"].copy()
            
            # Extract the numeric values in the () of 'V-REGION Nb of nucleotides'
            stats_IgE_mut_1.loc[:, 'V-REGION Nb of nucleotides-ex'] = stats_IgE_mut_1['V-REGION Nb of nucleotides'].str.replace(
            r'.*\((\d+)\)', r'\1', regex=True
            )
            # Extract the numeric values in the () of ''V-REGION Nb of mutations'
            stats_IgE_mut_1.loc[:, 'V-REGION Nb of mutations-ex'] = stats_IgE_mut_1['V-REGION Nb of mutations'].str.replace(
            r'.*\((\d+)\)', r'\1', regex=True
            )
            # Extract the numeric values in the () of ''CDR3-IMGT Nb of nonsilent mutations'
            stats_IgE_mut_1.loc[:, 'CDR3-IMGT Nb of nonsilent mutations-ex'] = stats_IgE_mut_1['CDR3-IMGT Nb of nonsilent mutations'].str.replace(
            r'.*\((\d+)\)', r'\1', regex=True
            )
            
            # Calculate the percentage of non-silent mutation in the CDR2 and CDR3 region
            stats_IgE_mut_2 = stats_IgE_mut_1.copy()
            
            # Ensure columns are numeric, coercing invalid parsing to NaN
            stats_IgE_mut_2['V-REGION Nb of nucleotides-ex'] = pd.to_numeric(stats_IgE_mut_2['V-REGION Nb of nucleotides-ex'], errors='coerce')
            stats_IgE_mut_2['CDR2-IMGT Nb of nonsilent mutations'] = pd.to_numeric(stats_IgE_mut_2['CDR2-IMGT Nb of nonsilent mutations'], errors='coerce')
            stats_IgE_mut_2['CDR3-IMGT Nb of nonsilent mutations-ex'] = pd.to_numeric(stats_IgE_mut_2['CDR3-IMGT Nb of nonsilent mutations-ex'], errors='coerce')
            
            # Calculate the percentage in the CDR2
            stats_IgE_mut_2['percentage of CDR2-IMGT Nb of nonsilent mutations'] = (
            stats_IgE_mut_2['CDR2-IMGT Nb of nonsilent mutations'] /
            stats_IgE_mut_2['V-REGION Nb of nucleotides-ex']
            ) * 100
            
            # Calculate the percentage in the CDR3
            stats_IgE_mut_2['percentage of CDR3-IMGT Nb of nonsilent mutations'] = (
            stats_IgE_mut_2['CDR3-IMGT Nb of nonsilent mutations-ex'] /
            stats_IgE_mut_2['V-REGION Nb of nucleotides-ex']
            ) * 100
            
            
            # Compare seqid of the IgE mutation table to the seqid of IgE-bias, IgG1-bias and unique IgE
            
            stats_IgE_mut_3 = stats_IgE_mut_2.copy()
            stats_IgE_mut_3['Match_in_IgE_biased_seqid'] = stats_IgE_mut_3['Sequence ID'].apply(lambda x: x in set(IgE_biased_clones_VDJ_compare['seqid']))
            stats_IgE_mut_3['Match_in_IgG1_biased_seqid'] = stats_IgE_mut_3['Sequence ID'].apply(lambda x: x in set(IgG1_biased_clones_VDJ_compare['seqid']))
            stats_IgE_mut_3['Match_in_unique_IgE_seqid'] = stats_IgE_mut_3['Sequence ID'].apply(lambda x: x in set(Unique_IgE['seqid']))
            stats_IgE_mut_3['Match_in_Total_IgE_seqid'] = stats_IgE_mut_3['Sequence ID'].apply(lambda x: x in set(stats_IgE_productive_VDJ['seqid']))
            
            # Filter out mutated and unmutated IgE-bias clones. mutated and unmutated IgG1-bias clones and mutated and unmutated unique IgE clones
            
            # Convert 'V-REGION Nb of mutations-ex' to numeric, coercing errors to NaN
            stats_IgE_mut_3['V-REGION Nb of mutations-ex'] = pd.to_numeric(stats_IgE_mut_3['V-REGION Nb of mutations-ex'], errors='coerce')
            
            # Ensure 'Match_in_IgE_biased_seqid' is a string and remove leading/trailing spaces
            stats_IgE_mut_3['Match_in_IgE_biased_seqid'] = stats_IgE_mut_3['Match_in_IgE_biased_seqid'].astype(str).str.strip()
            stats_IgE_mut_3['Match_in_IgG1_biased_seqid'] = stats_IgE_mut_3['Match_in_IgG1_biased_seqid'].astype(str).str.strip()
            stats_IgE_mut_3['Match_in_unique_IgE_seqid'] = stats_IgE_mut_3['Match_in_unique_IgE_seqid'].astype(str).str.strip()
            stats_IgE_mut_3['Match_in_Total_IgE_seqid'] = stats_IgE_mut_3['Match_in_Total_IgE_seqid'].astype(str).str.strip()
            
            # Now filter the mutated IgE_bias, IgG1_bias,Total IgE and unique IgE
            stats_tab_mutated_IgE_bias = stats_IgE_mut_3[(stats_IgE_mut_3['Match_in_IgE_biased_seqid'] == 'True') & (stats_IgE_mut_3['V-REGION Nb of mutations-ex'] > 0)]
            stats_tab_mutated_IgG1_bias = stats_IgE_mut_3[(stats_IgE_mut_3['Match_in_IgG1_biased_seqid'] == 'True') & (stats_IgE_mut_3['V-REGION Nb of mutations-ex'] > 0)]
            stats_tab_mutated_Unique_IgE = stats_IgE_mut_3[(stats_IgE_mut_3['Match_in_unique_IgE_seqid'] == 'True') & (stats_IgE_mut_3['V-REGION Nb of mutations-ex'] > 0)]
            stats_tab_mutated_Total_IgE = stats_IgE_mut_3[(stats_IgE_mut_3['Match_in_Total_IgE_seqid'] == 'True') & (stats_IgE_mut_3['V-REGION Nb of mutations-ex'] > 0)]
            
            # Now filter the unmutated IgE_bias, IgG1_bias,Total IgE and unique IgE
            stats_tab_unmutated_IgE_bias = stats_IgE_mut_3[(stats_IgE_mut_3['Match_in_IgE_biased_seqid'] == 'True') & (stats_IgE_mut_3['V-REGION Nb of mutations-ex'] == 0)]
            stats_tab_unmutated_IgG1_bias = stats_IgE_mut_3[(stats_IgE_mut_3['Match_in_IgG1_biased_seqid'] == 'True') & (stats_IgE_mut_3['V-REGION Nb of mutations-ex'] == 0)]
            stats_tab_unmutated_Unique_IgE = stats_IgE_mut_3[(stats_IgE_mut_3['Match_in_unique_IgE_seqid'] == 'True') & (stats_IgE_mut_3['V-REGION Nb of mutations-ex'] == 0)]
            stats_tab_unmutated_Total_IgE = stats_IgE_mut_3[(stats_IgE_mut_3['Match_in_Total_IgE_seqid'] == 'True') & (stats_IgE_mut_3['V-REGION Nb of mutations-ex'] == 0)]
            
            
            #Evaluate the mean of the percentage of CDR2 and CDR3 non-silent mutation of mutated IgE_bias, IgG1_bias,Total IgE and unique IgE
            mean_CDR2_stats_tab_mutated_IgE_bias = stats_tab_mutated_IgE_bias['percentage of CDR2-IMGT Nb of nonsilent mutations'].mean()
            print(f'Mutated IgE_biased_average percentage of CDR2 non silent mutation: {mean_CDR2_stats_tab_mutated_IgE_bias}')
            
            mean_CDR3_stats_tab_mutated_IgE_bias = stats_tab_mutated_IgE_bias['percentage of CDR3-IMGT Nb of nonsilent mutations'].mean()
            print(f'Mutated IgE_biased_average percentage of CDR3 non silent mutation: {mean_CDR3_stats_tab_mutated_IgE_bias}')
            
            mean_CDR2_stats_tab_mutated_IgG1_bias = stats_tab_mutated_IgG1_bias['percentage of CDR2-IMGT Nb of nonsilent mutations'].mean()
            print(f'Mutated IgG1_biased_average percentage of CDR2 non silent mutation: {mean_CDR2_stats_tab_mutated_IgG1_bias}')
            
            mean_CDR3_stats_tab_mutated_IgG1_bias = stats_tab_mutated_IgG1_bias['percentage of CDR3-IMGT Nb of nonsilent mutations'].mean()
            print(f'Mutated IgG1_biased_average percentage of CDR3 non silent mutation: {mean_CDR3_stats_tab_mutated_IgG1_bias}')
            
            mean_CDR2_stats_tab_mutated_Unique_IgE = stats_tab_mutated_Unique_IgE['percentage of CDR2-IMGT Nb of nonsilent mutations'].mean()
            print(f'Mutated unique_IgE_average percentage of CDR2 non silent mutation: {mean_CDR2_stats_tab_mutated_Unique_IgE}')
            
            mean_CDR3_stats_tab_mutated_Unique_IgE = stats_tab_mutated_Unique_IgE['percentage of CDR3-IMGT Nb of nonsilent mutations'].mean()
            print(f'Mutated unique_IgE__average percentage of CDR3 non silent mutation: {mean_CDR3_stats_tab_mutated_Unique_IgE}')
            
            #Evaluate the mean of the percentage of CDR2 and CDR3 non-silent mutation of unmutated IgE_bias, IgG1_bias,Total IgE and unique IgE
            if stats_tab_unmutated_IgE_bias.empty or stats_tab_unmutated_IgE_bias['percentage of CDR2-IMGT Nb of nonsilent mutations'].isna().all():
                mean_CDR2_stats_tab_unmutated_IgE_bias = 0
            else:
                mean_CDR2_stats_tab_unmutated_IgE_bias = stats_tab_unmutated_IgE_bias['percentage of CDR2-IMGT Nb of nonsilent mutations'].mean()
            print(f'Unmutated IgE_biased average percentage of CDR2 non silent mutation: {mean_CDR2_stats_tab_unmutated_IgE_bias}')
            

            if stats_tab_unmutated_IgE_bias.empty or stats_tab_unmutated_IgE_bias['percentage of CDR3-IMGT Nb of nonsilent mutations'].isna().all():
                mean_CDR3_stats_tab_unmutated_IgE_bias = 0
            else:
                mean_CDR3_stats_tab_unmutated_IgE_bias = stats_tab_unmutated_IgE_bias['percentage of CDR3-IMGT Nb of nonsilent mutations'].mean()
            print(f'Unmutated IgE_biased average percentage of CDR3 non silent mutation: {mean_CDR3_stats_tab_unmutated_IgE_bias}')
            

            if stats_tab_unmutated_IgG1_bias.empty or stats_tab_unmutated_IgG1_bias['percentage of CDR2-IMGT Nb of nonsilent mutations'].isna().all():
                mean_CDR2_stats_tab_unmutated_IgG1_bias = 0
            else:
                mean_CDR2_stats_tab_unmutated_IgG1_bias = stats_tab_unmutated_IgG1_bias['percentage of CDR2-IMGT Nb of nonsilent mutations'].mean()
            print(f'Unmutated IgG1_biased average percentage of CDR2 non silent mutation: {mean_CDR2_stats_tab_unmutated_IgG1_bias}')
            
            
            if stats_tab_unmutated_IgG1_bias.empty or stats_tab_unmutated_IgG1_bias['percentage of CDR3-IMGT Nb of nonsilent mutations'].isna().all():
                mean_CDR3_stats_tab_unmutated_IgG1_bias = 0
            else:
                mean_CDR3_stats_tab_unmutated_IgG1_bias = stats_tab_unmutated_IgG1_bias['percentage of CDR3-IMGT Nb of nonsilent mutations'].mean()
            print(f'Unmutated IgG1_biased average percentage of CDR3 non silent mutation: {mean_CDR3_stats_tab_unmutated_IgG1_bias}')
            

            if stats_tab_unmutated_Unique_IgE.empty or stats_tab_unmutated_Unique_IgE['percentage of CDR2-IMGT Nb of nonsilent mutations'].isna().all():
                mean_CDR2_stats_tab_unmutated_Unique_IgE = 0
            else:
                mean_CDR2_stats_tab_unmutated_Unique_IgE = stats_tab_unmutated_Unique_IgE['percentage of CDR2-IMGT Nb of nonsilent mutations'].mean()
            print(f'Unmutated Unique_IgE average percentage of CDR2 non silent mutation: {mean_CDR2_stats_tab_unmutated_Unique_IgE}')
            
            
            if stats_tab_unmutated_Unique_IgE.empty or stats_tab_unmutated_Unique_IgE['percentage of CDR3-IMGT Nb of nonsilent mutations'].isna().all():
                mean_CDR3_stats_tab_unmutated_Unique_IgE = 0
            else:
                mean_CDR3_stats_tab_unmutated_Unique_IgE = stats_tab_unmutated_Unique_IgE['percentage of CDR3-IMGT Nb of nonsilent mutations'].mean()
            print(f'Unmutated Unique_IgE average percentage of CDR3 non silent mutation: {mean_CDR3_stats_tab_unmutated_Unique_IgE}')
            
            #Barplot of the mean percentage CDR2 and CDR3 mutation
            # Data
            data = {
            'Subset': [
            'Mutated IgE-bias', 'Mutated IgE-bias',
            'Mutated IgG1-bias', 'Mutated IgG1-bias',
            'Mutated Unique IgE', 'Mutated Unique IgE',
            'Unmutated IgE-bias', 'Unmutated IgE-bias',
            'Unmutated IgG1-bias', 'Unmutated IgG1-bias',
            'Unmutated Unique IgE', 'Unmutated Unique IgE'
            ],
            'Region': ['CDR2', 'CDR3'] * 6,
            'Percentage': [
            mean_CDR2_stats_tab_mutated_IgE_bias, mean_CDR3_stats_tab_mutated_IgE_bias,
            mean_CDR2_stats_tab_mutated_IgG1_bias, mean_CDR3_stats_tab_mutated_IgG1_bias,
            mean_CDR2_stats_tab_mutated_Unique_IgE, mean_CDR3_stats_tab_mutated_Unique_IgE] +
            [mean_CDR2_stats_tab_unmutated_IgE_bias,mean_CDR3_stats_tab_unmutated_IgE_bias,
            mean_CDR2_stats_tab_unmutated_IgG1_bias, mean_CDR3_stats_tab_unmutated_IgG1_bias,
            mean_CDR2_stats_tab_unmutated_Unique_IgE, mean_CDR3_stats_tab_unmutated_Unique_IgE]
            
            }
            
            
            # Create DataFrame
            plot2mut = pd.DataFrame(data)
            
            # Plot
            fig4 = plt.figure(figsize=(10, 6))
            sns.barplot(data=plot2mut, x='Region', y='Percentage', hue='Subset', palette='deep')
            
            # Customize
            plt.title('Non-Silent Mutations in the CDR Regions', fontsize=14)
            plt.xlabel('VH Region', fontsize=12)
            plt.ylabel('Mean percentage non silent mutation(%)', fontsize=12)
            plt.legend(title='Subset', fontsize=10)
            plt.ylim(0, 0.8)  # Adjust limit if necessary
            plt.tight_layout()
            
            downloads_folder = os.path.join(os.path.expanduser("~"), "Downloads")
            
            plots_folder = os.path.join(downloads_folder, "BCR_Analysis_Plots")
            if not os.path.exists(plots_folder):
                os.makedirs(plots_folder)  # Create the folder if it doesn't exist
                
            plot_path = os.path.join(plots_folder, "Mean percentage non-silent mutation in the CDR2 and CDR3.png")
            plt.savefig(plot_path, dpi=1200, bbox_inches='tight')
                        
            # Create a copy of the dataframe
            stats_IgE_productive_VDJ_seqid_compare = stats_IgE_productive_VDJ.copy()
            
            #Compare the seqid from Total IgE (containing the VDJ) to the seqid of mutated IgE_bias, IgG1_bias,Total IgE and unique IgE
            stats_IgE_productive_VDJ_seqid_compare['Match_in_mutated_IgE_biased_seqid'] = stats_IgE_productive_VDJ_seqid_compare['seqid'].apply(lambda x: x in set(stats_tab_mutated_IgE_bias['Sequence ID']))
            stats_IgE_productive_VDJ_seqid_compare['Match_in_mutated_IgG1_biased_seqid'] = stats_IgE_productive_VDJ_seqid_compare['seqid'].apply(lambda x: x in set(stats_tab_mutated_IgG1_bias['Sequence ID']))
            stats_IgE_productive_VDJ_seqid_compare['Match_in_mutated_unique_IgE_seqid'] = stats_IgE_productive_VDJ_seqid_compare['seqid'].apply(lambda x: x in set(stats_tab_mutated_Unique_IgE['Sequence ID']))
            stats_IgE_productive_VDJ_seqid_compare['Match_in_mutated_Total_IgE_seqid'] = stats_IgE_productive_VDJ_seqid_compare['seqid'].apply(lambda x: x in set(stats_tab_mutated_Total_IgE['Sequence ID']))
            
            #Compare the seqid from Total IgE (containing the VDJ) to the seqid of unmutated IgE_bias, IgG1_bias,Total IgE and unique IgE
            stats_IgE_productive_VDJ_seqid_compare['Match_in_unmutated_IgE_biased_seqid'] = stats_IgE_productive_VDJ_seqid_compare['seqid'].apply(lambda x: x in set(stats_tab_unmutated_IgE_bias['Sequence ID']))
            stats_IgE_productive_VDJ_seqid_compare['Match_in_unmutated_IgG1_biased_seqid'] = stats_IgE_productive_VDJ_seqid_compare['seqid'].apply(lambda x: x in set(stats_tab_unmutated_IgG1_bias['Sequence ID']))
            stats_IgE_productive_VDJ_seqid_compare['Match_in_unmutated_unique_IgE_seqid'] = stats_IgE_productive_VDJ_seqid_compare['seqid'].apply(lambda x: x in set(stats_tab_unmutated_Unique_IgE['Sequence ID']))
            stats_IgE_productive_VDJ_seqid_compare['Match_in_unmutated_Total_IgE_seqid'] = stats_IgE_productive_VDJ_seqid_compare['seqid'].apply(lambda x: x in set(stats_tab_unmutated_Total_IgE['Sequence ID']))
            
            
            # Ensure 'Match_in_mutated IgE_bias, IgG1_bias,Total IgE and unique IgE seqid' is a string and remove leading/trailing spaces
            stats_IgE_productive_VDJ_seqid_compare['Match_in_mutated_IgE_biased_seqid'] = stats_IgE_productive_VDJ_seqid_compare['Match_in_mutated_IgE_biased_seqid'].astype(str).str.strip()
            stats_IgE_productive_VDJ_seqid_compare['Match_in_mutated_IgG1_biased_seqid'] = stats_IgE_productive_VDJ_seqid_compare['Match_in_mutated_IgG1_biased_seqid'].astype(str).str.strip()
            stats_IgE_productive_VDJ_seqid_compare['Match_in_mutated_unique_IgE_seqid'] = stats_IgE_productive_VDJ_seqid_compare['Match_in_mutated_unique_IgE_seqid'].astype(str).str.strip()
            stats_IgE_productive_VDJ_seqid_compare['Match_in_mutated_Total_IgE_seqid'] = stats_IgE_productive_VDJ_seqid_compare['Match_in_mutated_Total_IgE_seqid'].astype(str).str.strip()
            
            # Ensure 'Match_in_unmutated IgE_bias, IgG1_bias,Total IgE and unique IgE seqid' is a string and remove leading/trailing spaces
            stats_IgE_productive_VDJ_seqid_compare['Match_in_unmutated_IgE_biased_seqid'] = stats_IgE_productive_VDJ_seqid_compare['Match_in_unmutated_IgE_biased_seqid'].astype(str).str.strip()
            stats_IgE_productive_VDJ_seqid_compare['Match_in_unmutated_IgG1_biased_seqid'] = stats_IgE_productive_VDJ_seqid_compare['Match_in_unmutated_IgG1_biased_seqid'].astype(str).str.strip()
            stats_IgE_productive_VDJ_seqid_compare['Match_in_unmutated_unique_IgE_seqid'] = stats_IgE_productive_VDJ_seqid_compare['Match_in_unmutated_unique_IgE_seqid'].astype(str).str.strip()
            stats_IgE_productive_VDJ_seqid_compare['Match_in_unmutated_Total_IgE_seqid'] = stats_IgE_productive_VDJ_seqid_compare['Match_in_unmutated_Total_IgE_seqid'].astype(str).str.strip()
            
            
            # Now filter the mutated IgE_bias, IgG1_bias,Total IgE and unique IgE clones
            stats_tab_mutated_IgE_bias_VDJ = stats_IgE_productive_VDJ_seqid_compare[(stats_IgE_productive_VDJ_seqid_compare['Match_in_mutated_IgE_biased_seqid'] == 'True')]
            stats_tab_mutated_IgG1_bias_VDJ = stats_IgE_productive_VDJ_seqid_compare[(stats_IgE_productive_VDJ_seqid_compare['Match_in_mutated_IgG1_biased_seqid'] == 'True')]
            stats_tab_mutated_unique_IgE_VDJ = stats_IgE_productive_VDJ_seqid_compare[(stats_IgE_productive_VDJ_seqid_compare['Match_in_mutated_unique_IgE_seqid'] == 'True')]
            stats_tab_mutated_Total_IgE_VDJ = stats_IgE_productive_VDJ_seqid_compare[(stats_IgE_productive_VDJ_seqid_compare['Match_in_mutated_Total_IgE_seqid'] == 'True')]
            
            # Now filter the unmutated IgE_bias, IgG1_bias,Total IgE and unique IgE clones
            stats_tab_unmutated_IgE_bias_VDJ = stats_IgE_productive_VDJ_seqid_compare[(stats_IgE_productive_VDJ_seqid_compare['Match_in_unmutated_IgE_biased_seqid'] == 'True')]
            stats_tab_unmutated_IgG1_bias_VDJ = stats_IgE_productive_VDJ_seqid_compare[(stats_IgE_productive_VDJ_seqid_compare['Match_in_unmutated_IgG1_biased_seqid'] == 'True')]
            stats_tab_unmutated_unique_IgE_VDJ = stats_IgE_productive_VDJ_seqid_compare[(stats_IgE_productive_VDJ_seqid_compare['Match_in_unmutated_unique_IgE_seqid'] == 'True')]
            stats_tab_unmutated_Total_IgE_VDJ = stats_IgE_productive_VDJ_seqid_compare[(stats_IgE_productive_VDJ_seqid_compare['Match_in_unmutated_Total_IgE_seqid'] == 'True')]
            
            
            
            # Number and percentage of mutated IgE_bias, IgG1_bias,Total IgE and unique IgE clones (Unique VDJ), average copy number of all IgE subsets
            unique_strings_mutated_IgE_bias_VDJ = stats_tab_mutated_IgE_bias_VDJ['IgE_VDJ'].nunique(dropna=True)
            print(f'Number of mutated IgE-biased clones: {unique_strings_mutated_IgE_bias_VDJ}')
            
            percentage_mutated_IgE_bias_clones = (unique_strings_mutated_IgE_bias_VDJ / unique_strings_IgE_productive_VDJ) * 100
            print(f'percentage of mutated IgE-biased clones: {percentage_mutated_IgE_bias_clones}')
            
            mean_copynumber_mutated_IgE_bias_VDJ = stats_tab_mutated_IgE_bias_VDJ['onecopy'].mean()
            print(f'mean copy number of mutated IgE-biased clones: {mean_copynumber_mutated_IgE_bias_VDJ}')
            
            
            unique_strings_mutated_IgG1_bias_VDJ = stats_tab_mutated_IgG1_bias_VDJ['IgE_VDJ'].nunique(dropna=True)
            print(f'Number of mutated IgG1-biased clones: {unique_strings_mutated_IgG1_bias_VDJ}')
            
            percentage_mutated_IgG1_bias_clones = (unique_strings_mutated_IgG1_bias_VDJ / unique_strings_IgE_productive_VDJ) * 100
            print(f'percentage of mutated IgG1-biased clones: {percentage_mutated_IgG1_bias_clones}')
            
            mean_copynumber_mutated_IgG1_bias_VDJ = stats_tab_mutated_IgG1_bias_VDJ['onecopy'].mean()
            print(f'mean copy number of mutated IgG1-biased clones: {mean_copynumber_mutated_IgG1_bias_VDJ}')
            
            
            unique_strings_mutated_unique_IgE_VDJ = stats_tab_mutated_unique_IgE_VDJ['IgE_VDJ'].nunique(dropna=True)
            print(f'Number of mutated unique IgE clones: {unique_strings_mutated_unique_IgE_VDJ}')
            
            percentage_mutated_unique_IgE_clones = (unique_strings_mutated_unique_IgE_VDJ / unique_strings_IgE_productive_VDJ) * 100
            print(f'percentage of mutated unique IgE clones: {percentage_mutated_unique_IgE_clones}')
            
            mean_copynumber_mutated_unique_IgE_VDJ = stats_tab_mutated_unique_IgE_VDJ['onecopy'].mean()
            print(f'mean copy number of mutated unique IgE clones: {mean_copynumber_mutated_unique_IgE_VDJ}')
            
            
            unique_strings_unmutated_Total_IgE_VDJ = stats_tab_unmutated_Total_IgE_VDJ['IgE_VDJ'].nunique(dropna=True)
            
            unique_strings_mutated_Total_IgE_VDJ = stats_tab_mutated_Total_IgE_VDJ['IgE_VDJ'].nunique(dropna=True)
            print(f'Number of mutated Total IgE clones: {unique_strings_mutated_Total_IgE_VDJ}')
            
            try:
                percentage_mutated_Total_IgE_clones = (unique_strings_mutated_Total_IgE_VDJ / (unique_strings_mutated_Total_IgE_VDJ + unique_strings_unmutated_Total_IgE_VDJ)) * 100
                print(f'Percentage of mutated Total IgE clones: {percentage_mutated_Total_IgE_clones}')
            except ZeroDivisionError:
                print('Percentage of mutated Total IgE clones: 0% (No clones present)')
            
            # Define the function to print the result
            def your_script_function(mean_copynumber_mutated_Total_IgE_VDJ):
                print(f'mean copy number of mutated Total IgE clones: {mean_copynumber_mutated_Total_IgE_VDJ}')
            
            # Calculate the value outside the function
            mean_copynumber_mutated_Total_IgE_VDJ = stats_tab_mutated_Total_IgE_VDJ['onecopy'].mean()
            
            # Call the function with the calculated value
            your_script_function(mean_copynumber_mutated_Total_IgE_VDJ)
                
            # Number and percentage of mutated IgE_bias, IgG1_bias,Total IgE and unique IgE clones (Unique VDJ), average copy number of all IgE subsets
            
            unique_strings_unmutated_IgE_bias_VDJ = stats_tab_unmutated_IgE_bias_VDJ['IgE_VDJ'].nunique(dropna=True)
            print(f'Number of unmutated IgE-biased clones: {unique_strings_unmutated_IgE_bias_VDJ}')
            
            percentage_unmutated_IgE_bias_clones = (unique_strings_unmutated_IgE_bias_VDJ / unique_strings_IgE_productive_VDJ) * 100
            print(f'percentage of unmutated IgE-biased clones: {percentage_unmutated_IgE_bias_clones}')
            
            mean_copynumber_unmutated_IgE_bias_VDJ = stats_tab_unmutated_IgE_bias_VDJ['onecopy'].mean()
            
            mean_copynumber_unmutated_IgE_bias_VDJ = 0 if pd.isna(mean_copynumber_unmutated_IgE_bias_VDJ) else mean_copynumber_unmutated_IgE_bias_VDJ
            
            print(f'mean copy number of unmutated IgE-biased clones: {mean_copynumber_unmutated_IgE_bias_VDJ}')
            
            
            unique_strings_unmutated_IgG1_bias_VDJ = stats_tab_unmutated_IgG1_bias_VDJ['IgE_VDJ'].nunique(dropna=True)
            print(f'Number of unmutated IgG1-biased clones: {unique_strings_unmutated_IgG1_bias_VDJ}')
            
            percentage_unmutated_IgG1_bias_clones = (unique_strings_unmutated_IgG1_bias_VDJ / unique_strings_IgE_productive_VDJ) * 100
            print(f'percentage of unmutated IgG1-biased clones: {percentage_unmutated_IgG1_bias_clones}')
            
            mean_copynumber_unmutated_IgG1_bias_VDJ = stats_tab_unmutated_IgG1_bias_VDJ['onecopy'].mean()
            
            mean_copynumber_unmutated_IgG1_bias_VDJ = 0 if pd.isna(mean_copynumber_unmutated_IgG1_bias_VDJ) else mean_copynumber_unmutated_IgG1_bias_VDJ
            
            print(f'mean copy number of unmutated IgG1-biased clones: {mean_copynumber_unmutated_IgG1_bias_VDJ}')
            
            
            unique_strings_unmutated_unique_IgE_VDJ = stats_tab_unmutated_unique_IgE_VDJ['IgE_VDJ'].nunique(dropna=True)
            print(f'Number of unmutated unique IgE clones: {unique_strings_unmutated_unique_IgE_VDJ}')
            
            percentage_unmutated_unique_IgE_clones = (unique_strings_unmutated_unique_IgE_VDJ / unique_strings_IgE_productive_VDJ) * 100
            print(f'percentage of unmutated unique IgE clones: {percentage_unmutated_unique_IgE_clones}')
            
            mean_copynumber_unmutated_unique_IgE_VDJ = stats_tab_unmutated_unique_IgE_VDJ['onecopy'].mean()
            
            mean_copynumber_unmutated_unique_IgE_VDJ = 0 if pd.isna(mean_copynumber_unmutated_unique_IgE_VDJ) else mean_copynumber_unmutated_unique_IgE_VDJ
            
            print(f'mean copy number of unmutated unique IgE clones: {mean_copynumber_unmutated_unique_IgE_VDJ}')
            
            
            print(f'Number of unmutated Total IgE clones: {unique_strings_unmutated_Total_IgE_VDJ}')
            
            total_unique_strings = unique_strings_mutated_Total_IgE_VDJ + unique_strings_unmutated_Total_IgE_VDJ
            if total_unique_strings > 0:
                percentage_unmutated_Total_IgE_clones = (unique_strings_unmutated_Total_IgE_VDJ / total_unique_strings) * 100
                print(f'Percentage of unmutated Total IgE clones: {percentage_unmutated_Total_IgE_clones}')
                print('Percentage of unmutated Total IgE clones: 0 (0)')
            
            
            
            mean_copynumber_unmutated_Total_IgE_VDJ = stats_tab_unmutated_Total_IgE_VDJ['onecopy'].mean()
            mean_copynumber_unmutated_Total_IgE_VDJ = 0 if pd.isna(mean_copynumber_unmutated_Total_IgE_VDJ) else mean_copynumber_unmutated_Total_IgE_VDJ
            
            print(f'mean copy number of unmutated Total IgE clones: {mean_copynumber_unmutated_Total_IgE_VDJ}')
            
            
            # Cleanup IgG1 mutation table, remove Non-productive clones and cleanup the number of mutations

            # Check if either 'V.DOMAIN.Functionality' or 'V-DOMAIN.Functionality' exists in the DataFrame
            if 'V.DOMAIN.Functionality' in stats_IgG1_mut.columns:
                column_name = 'V.DOMAIN.Functionality'
            elif 'V-DOMAIN Functionality' in stats_IgG1_mut.columns:
                column_name = 'V-DOMAIN Functionality'
            else:
                raise KeyError("Neither 'V.DOMAIN.Functionality' nor 'V-DOMAIN Functionality' column is present in the DataFrame.")
            
            # Filter for productive sequences
            stats_IgG1_mut_1 = stats_IgG1_mut[stats_IgG1_mut[column_name] == "productive"].copy()
            
            # Sort the filtered DataFrame by the selected column
            stats_IgG1_mut_1_sorted = stats_IgG1_mut_1.sort_values(by=column_name)
            
            
            # Check for the correct column name
            if 'V.REGION.Nb.of.nucleotides' in stats_IgG1_mut_1.columns:
                column_name = 'V.REGION.Nb.of.nucleotides'
            elif 'V-REGION Nb of nucleotides' in stats_IgG1_mut_1.columns:
                column_name = 'V-REGION Nb of nucleotides'
            else:
                raise KeyError("Neither 'V.REGION.Nb.of.nucleotides' nor 'V-REGION Nb of nucleotides' column is present in the DataFrame.")
            
            # Apply str.replace to extract the number within parentheses
            stats_IgG1_mut_1['V-REGION Nb of nucleotides-ex'] = stats_IgG1_mut_1[column_name].str.replace(
                r'.*\((\d+)\)', r'\1', regex=True
            )

            
            # Handle 'V.REGION.Nb.of.mutations' or 'V-REGION Nb of mutations'
            if 'V.REGION.Nb.of.mutations' in stats_IgG1_mut_1.columns:
                column_mutations = 'V.REGION.Nb.of.mutations'
            elif 'V-REGION Nb of mutations' in stats_IgG1_mut_1.columns:
                column_mutations = 'V-REGION Nb of mutations'
            else:
                raise KeyError("Neither 'V.REGION.Nb.of.mutations' nor 'V-REGION Nb of mutations' column is present in the DataFrame.")
            
            stats_IgG1_mut_1['V-REGION Nb of mutations-ex'] = stats_IgG1_mut_1[column_mutations].str.replace(
                r'.*\((\d+)\)', r'\1', regex=True
            )


            # Handle 'CDR3.IMGT.Nb.of.nonsilent.mutations' or 'CDR3-IMGT Nb of nonsilent mutations'
            if 'CDR3.IMGT.Nb.of.nonsilent.mutations' in stats_IgG1_mut_1.columns:
                column_cdr3 = 'CDR3.IMGT.Nb.of.nonsilent.mutations'
            elif 'CDR3-IMGT Nb of nonsilent mutations' in stats_IgG1_mut_1.columns:
                column_cdr3 = 'CDR3-IMGT Nb of nonsilent mutations'
            else:
                raise KeyError("Neither 'CDR3.IMGT.Nb.of.nonsilent.mutations' nor 'CDR3-IMGT Nb of nonsilent mutations' column is present in the DataFrame.")
            
            # Create the new column with a consistent name
            stats_IgG1_mut_1['CDR3-IMGT Nb of nonsilent mutations-ex'] = stats_IgG1_mut_1[column_cdr3].str.replace(
                r'.*\((\d+)\)', r'\1', regex=True
            )
            
            
            # Ensure stats_IgG1_mut_1 is a Pandas DataFrame
            stats_IgG1_mut_3 = stats_IgG1_mut_1.copy()  # Convert Dask to Pandas DataFrame
            
            # Add a column to indicate if 'Sequence.ID' is in the productive IgG1 VDJ
            # Check for the correct column name
            if 'Sequence.ID' in stats_IgG1_mut_3.columns:
                column_seqid = 'Sequence.ID'
            elif 'Sequence ID' in stats_IgG1_mut_3.columns:
                column_seqid = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present in the DataFrame.")
            
            # Perform analysis using the detected column
            stats_IgG1_mut_3['Match_in_Total_IgG1_seqid'] = stats_IgG1_mut_3[column_seqid].isin(stats_IgG1_productive_VDJ['seqid'])

            
            # Convert 'V-REGION Nb of mutations-ex' to numeric, handling errors by coercing invalid entries
            stats_IgG1_mut_3['V-REGION Nb of mutations-ex'] = pd.to_numeric(stats_IgG1_mut_3['V-REGION Nb of mutations-ex'], errors='coerce')
            
            # Filter mutated and unmutated Total IgG1 directly using boolean masks
            stats_tab_mutated_Total_IgG1 = stats_IgG1_mut_3[
            (stats_IgG1_mut_3['Match_in_Total_IgG1_seqid']) & (stats_IgG1_mut_3['V-REGION Nb of mutations-ex'] > 0)
            ]
            
            stats_tab_unmutated_Total_IgG1 = stats_IgG1_mut_3[
            (stats_IgG1_mut_3['Match_in_Total_IgG1_seqid']) & (stats_IgG1_mut_3['V-REGION Nb of mutations-ex'] == 0)
            ]
            
            # Create sets of seqids for faster lookup

            # Check for the correct column name in stats_tab_mutated_Total_IgG1
            if 'Sequence.ID' in stats_tab_mutated_Total_IgG1.columns:
                column_seqid = 'Sequence.ID'
            elif 'Sequence ID' in stats_tab_mutated_Total_IgG1.columns:
                column_seqid = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present in the DataFrame.")
            
            # Create a set of unique sequence IDs from the detected column
            mutated_seqids = set(stats_tab_mutated_Total_IgG1[column_seqid])

            
            # Check for the correct column name in stats_tab_unmutated_Total_IgG1
            if 'Sequence.ID' in stats_tab_unmutated_Total_IgG1.columns:
                column_seqid = 'Sequence.ID'
            elif 'Sequence ID' in stats_tab_unmutated_Total_IgG1.columns:
                column_seqid = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present in the DataFrame.")
            
            # Create a set of unique unmutated sequence IDs from the detected column
            unmutated_seqids = set(stats_tab_unmutated_Total_IgG1[column_seqid])


            # Create a copy of the dataframe for comparison
            stats_IgG1_productive_VDJ_seqid_compare = stats_IgG1_productive_VDJ.copy()
            
            # Update the comparison DataFrame with boolean masks
            stats_IgG1_productive_VDJ_seqid_compare['Match_in_mutated_Total_IgG1_seqid'] = stats_IgG1_productive_VDJ_seqid_compare['seqid'].isin(mutated_seqids)
            stats_IgG1_productive_VDJ_seqid_compare['Match_in_unmutated_Total_IgG1_seqid'] = stats_IgG1_productive_VDJ_seqid_compare['seqid'].isin(unmutated_seqids)
            
            
            # Ensure 'Match_in_mutated Total IgG1 seqid' is a string and remove leading/trailing spaces
            stats_IgG1_productive_VDJ_seqid_compare['Match_in_mutated_Total_IgG1_seqid'] = stats_IgG1_productive_VDJ_seqid_compare['Match_in_mutated_Total_IgG1_seqid'].astype(str).str.strip()
            
            # Ensure 'Match_in_unmutated Total IgG1 seqid' is a string and remove leading/trailing spaces
            stats_IgG1_productive_VDJ_seqid_compare['Match_in_unmutated_Total_IgG1_seqid'] = stats_IgG1_productive_VDJ_seqid_compare['Match_in_unmutated_Total_IgG1_seqid'].astype(str).str.strip()
            
            # Now filter the mutated Total IgG1 clones
            stats_tab_mutated_Total_IgG1_VDJ = stats_IgG1_productive_VDJ_seqid_compare[(stats_IgG1_productive_VDJ_seqid_compare['Match_in_mutated_Total_IgG1_seqid'] == 'True')]
            
            # Now filter the unmutated Total IgG1 clones
            stats_tab_unmutated_Total_IgG1_VDJ = stats_IgG1_productive_VDJ_seqid_compare[(stats_IgG1_productive_VDJ_seqid_compare['Match_in_unmutated_Total_IgG1_seqid'] == 'True')]
            
            
            # Number and percentage of mutated Total IgG1 clones, average copy number of mutated Total IgG1 clones
            
            # Calculate the number of unique strings in mutated Total IgG1 VDJ
            unique_strings_mutated_Total_IgG1_VDJ = stats_tab_mutated_Total_IgG1_VDJ['IgG1_VDJ'].nunique(dropna=True)
            print(f'Number of mutated Total IgG1 clones: {unique_strings_mutated_Total_IgG1_VDJ}')
            
            # Calculate the number of unique strings in unmutated Total IgG1 VDJ
            unique_strings_unmutated_Total_IgG1_VDJ = stats_tab_unmutated_Total_IgG1_VDJ['IgG1_VDJ'].nunique(dropna=True)
            print(f'Number of unmutated Total IgG1 clones: {unique_strings_unmutated_Total_IgG1_VDJ}')
            
            # Calculate the percentage of mutated Total IgG1 clones
            percentage_mutated_Total_IgG1_clones = (
            unique_strings_mutated_Total_IgG1_VDJ /
            (unique_strings_mutated_Total_IgG1_VDJ + unique_strings_unmutated_Total_IgG1_VDJ)
            ) * 100
            print(f'Percentage of mutated Total IgG1 clones: {percentage_mutated_Total_IgG1_clones:.2f}%')
            
            # Mean copy number of mutated Total IgG1 clones
            mean_copynumber_mutated_Total_IgG1_VDJ = stats_tab_mutated_Total_IgG1_VDJ['onecopy'].mean()
            print(f'Mean copy number of mutated Total IgG1 clones: {mean_copynumber_mutated_Total_IgG1_VDJ}')
            
            # Calculate the percentage of unmutated Total IgG1 clones
            percentage_unmutated_Total_IgG1_clones = (
            unique_strings_unmutated_Total_IgG1_VDJ /
            (unique_strings_mutated_Total_IgG1_VDJ + unique_strings_unmutated_Total_IgG1_VDJ)
            ) * 100
            print(f'Percentage of unmutated Total IgG1 clones: {percentage_unmutated_Total_IgG1_clones:.2f}%')
            
            # Mean copy number of unmutated Total IgG1 clones
            mean_copynumber_unmutated_Total_IgG1_VDJ = stats_tab_unmutated_Total_IgG1_VDJ['onecopy'].mean()
            print(f'Mean copy number of unmutated Total IgG1 clones: {mean_copynumber_unmutated_Total_IgG1_VDJ}')
            
            
            #Plot a pie chat of mutated and unmutated IgE and IgG1
            
            # Data for the first pie chart: Percentage of mutated IgE and IgG1
            percentage_mutated_IgE = percentage_mutated_Total_IgE_clones
    
            percentage_unmutated_IgE = percentage_unmutated_Total_IgE_clones
            
            labels1 = ['Mutated IgE', 'Unmutated IgE']
            sizes1 = [percentage_mutated_IgE, percentage_unmutated_IgE]
            colors1 = ['red', 'blue']
            explode1 = (0.1, 0)  # Highlight the first slice (mutated IgE)
            
            # Data for the second pie chart: Percentage of mutated and unmutated IgG1
            percentage_mutated_Total_IgG1 = percentage_mutated_Total_IgG1_clones
            percentage_unmutated_Total_IgG1 = percentage_unmutated_Total_IgG1_clones
            
            labels2 = ['Mutated IgG1', 'Unmutated IgG1']
            sizes2 = [percentage_mutated_Total_IgG1, percentage_unmutated_Total_IgG1]
            colors2 = ['red', 'blue']
            explode2 = (0.1, 0)  # Highlight the first slice (mutated IgG1)
            
            # Create subplots for the pie charts
            fig5, axis = plt.subplots(1, 2, figsize=(12, 6))
            
            # Plot the first pie chart
            axis[0].pie(sizes1, explode=explode1, labels=labels1, colors=colors1, autopct='%1.1f%%', startangle=90)
            axis[0].set_title('Percentage of Mutated and Unmutated IgE')
            
            # Plot the second pie chart
            axis[1].pie(sizes2, explode=explode2, labels=labels2, colors=colors2, autopct='%1.1f%%', startangle=90)
            axis[1].set_title('Percentage of Mutated and Unmutated IgG1')
            
            # Adjust layout and display the plots
            plt.tight_layout()
            downloads_folder = os.path.join(os.path.expanduser("~"), "Downloads")
            
            plots_folder = os.path.join(downloads_folder, "BCR_Analysis_Plots")
            if not os.path.exists(plots_folder):
                os.makedirs(plots_folder)  # Create the folder if it doesn't exist
                
            plot_path = os.path.join(plots_folder, "Percentage of mutated, unmutated IgE and IgG1.png")
            plt.savefig(plot_path, dpi=1200, bbox_inches='tight')
            
            
            
            
            #Bar plot number of clones from different IgE subsets
            types = ['IgE','IgE', 'IgE','IgE','IgE','IgE','IgE', 'IgE','IgG1','IgG1']
            categories = ['Mutated IgE', 'Unmutated IgE', 'Mutated IgE-biased', 'Unmutated IgE-biased',' Mutated IgG1-biased',' Unmutated IgG1-biased', 'Mutated unique IgE', 'Unmutated unique IgE',' Mutated IgG1', 'Unmutated IgG1']
            counts = [unique_strings_mutated_Total_IgE_VDJ, unique_strings_unmutated_Total_IgE_VDJ, unique_strings_mutated_IgE_bias_VDJ, unique_strings_unmutated_IgE_bias_VDJ,unique_strings_mutated_IgG1_bias_VDJ, unique_strings_unmutated_IgG1_bias_VDJ,unique_strings_mutated_unique_IgE_VDJ,unique_strings_unmutated_unique_IgE_VDJ,unique_strings_mutated_Total_IgG1_VDJ,unique_strings_unmutated_Total_IgG1_VDJ]
            
            # Create a DataFrame
            Graph3 = pd.DataFrame({
            'Type': types,
            'Category': categories,
            'Count': counts
            })
            # Set Seaborn style for scientific plots
            sns.set(style="whitegrid")
            
            # Filter data for IgE and IgG1
            df_IgE = Graph3[Graph3['Type'] == 'IgE']
            df_IgG1 = Graph3[Graph3['Type'] == 'IgG1']
            
            # Calculate dynamic y-axis limits
            max_IgE = df_IgE['Count'].max()
            max_IgG1 = df_IgG1['Count'].max()
            
            ylim_IgE = max_IgE + 10  # Slightly above the max value
            ylim_IgG1 = max_IgG1 + 100  # Slightly above the max value
            
            # Generate dynamic color palette based on the number of unique categories
            unique_categories_IgE = df_IgE['Category'].nunique()
            unique_categories_IgG1 = df_IgG1['Category'].nunique()
            
            # Get n unique colors for each plot
            palette_IgE = sns.color_palette("Set2", unique_categories_IgE)
            palette_IgG1 = sns.color_palette("Set2", unique_categories_IgG1)
            
            # Create two subplots
            fig6, axes = plt.subplots(1, 2, figsize=(12, 6))
            
            # Plot for IgE with dynamically generated color palette
            sns.barplot(
            data=df_IgE,
            x='Category',
            y='Count',
            hue='Category',
            palette=palette_IgE,
            dodge=False,
            ax=axes[0]
            )
            axes[0].set_title('Graph A: Mutated and Unmutated IgE subsets', fontsize=14)
            axes[0].set_xlabel('', fontsize=12)
            axes[0].set_ylabel('Number of Clones', fontsize=12)
            axes[0].set_ylim(0, ylim_IgE)  # Set dynamic y-axis limit for IgE
            axes[0].legend([], [], frameon=False)  # Remove legend
            axes[0].tick_params(axis='x', rotation=90)
            for container in axes[0].containers:
                axes[0].bar_label(container, fmt='%d', label_type='edge', fontsize=10)
                # Plot for IgG1 with dynamically generated color palette
                sns.barplot(
                data=df_IgG1,
                x='Category',
                y='Count',
                hue='Category',
                palette=palette_IgG1,
                dodge=False,
                ax=axes[1]
                )
                axes[1].set_title('Graph B: Mutated and Unmutated IgG1', fontsize=14)
                axes[1].set_xlabel('', fontsize=12)
                axes[1].set_ylabel('')  # No ylabel for the second plot
                axes[1].set_ylim(0, ylim_IgG1)  # Set dynamic y-axis limit for IgG1
                axes[1].legend([], [], frameon=False)  # Remove legend
                for container in axes[1].containers:
                    axes[1].bar_label(container, fmt='%d', label_type='edge', fontsize=10)
                
                # Adjust layout
                plt.tight_layout()
            
                downloads_folder = os.path.join(os.path.expanduser("~"), "Downloads")
            
            plots_folder = os.path.join(downloads_folder, "BCR_Analysis_Plots")
            if not os.path.exists(plots_folder):
                os.makedirs(plots_folder)  # Create the folder if it doesn't exist
                
            plot_path = os.path.join(plots_folder, "Number of mutated,unmutated IgE-bias,IgG1-bias and unique IgE and IgG1.png")
            plt.savefig(plot_path, dpi=1200, bbox_inches='tight')
            
            
            
            # Plot copy number distribution of mutated IgE and IgG1
            # Sort the counts in ascending order
            mutated_Total_IgE_copy_sort = stats_tab_mutated_Total_IgE_VDJ.sort_values(by='onecopy', ascending=True).reset_index(drop=True)
            mutated_Total_IgG1_copy_sort = stats_tab_mutated_Total_IgG1_VDJ.sort_values(by='onecopy', ascending=True).reset_index(drop=True)
            
            # Calculate the mean and median values for both dataframes
            mean1 = mutated_Total_IgE_copy_sort ['onecopy'].mean()
            mean2 = mutated_Total_IgG1_copy_sort ['onecopy'].mean()
            
            # Set Seaborn style
            sns.set(style="whitegrid")
            
            # Create the plot with two rows and two columns of subplots (original plots and zoomed-in plots)
            fig7, axes = plt.subplots(2, 2, figsize=(14, 12))
            
            # Plot for DataFrame 1 (Shared IgE) - Original
            sns.scatterplot(data=mutated_Total_IgE_copy_sort, x=mutated_Total_IgE_copy_sort.index, y='onecopy', color='blue', s=100, ax=axes[0, 0], label='Mutated IgE', edgecolor='w')
            axes[0, 0].axhline(mean1, color='orange', linestyle='--', label=f'IgE Mean: {mean1:.2f}')
            axes[0, 0].axhline(mean2, color='green', linestyle='--', label=f'IgG1 Mean: {mean2:.2f}')
            axes[0, 0].set_title('Clonal Size Distribution - Mutated IgE', fontsize=14)
            axes[0, 0].set_xlabel('Clonal Index', fontsize=12)
            axes[0, 0].set_ylabel('Clonal Size', fontsize=12)
            axes[0, 0].legend()
            
            # Plot for DataFrame 2 (Shared IgG1) - Original
            sns.scatterplot(data=mutated_Total_IgG1_copy_sort, x=mutated_Total_IgG1_copy_sort.index, y='onecopy', color='red', s=100, ax=axes[0, 1], label='Mutated IgG1', edgecolor='w')
            axes[0, 1].axhline(mean1, color='orange', linestyle='--', label=f'IgE Mean: {mean1:.2f}')
            axes[0, 1].axhline(mean2, color='green', linestyle='--', label=f'IgG1 Mean: {mean2:.2f}')
            axes[0, 1].set_title('Clonal Size Distribution - Mutated IgG1', fontsize=14)
            axes[0, 1].set_xlabel('Clonal Index', fontsize=12)
            axes[0, 1].set_ylabel('Clonal Size', fontsize=12)
            axes[0, 1].legend()
            
            # Define a range around the mean for zooming in (e.g., +-10 around the mean)
            zoom_range_1 = (mean2 - 20, mean2 + 20)
            zoom_range_2 = (mean2 - 20, mean2 + 20)
            
            # Filter DataFrames to show values around the mean
            Mutated_IgE_zoom = mutated_Total_IgE_copy_sort[(mutated_Total_IgE_copy_sort['onecopy'] >= zoom_range_1[0]) & (mutated_Total_IgE_copy_sort['onecopy'] <= zoom_range_1[1])]
            Mutated_IgG1_zoom = mutated_Total_IgG1_copy_sort[(mutated_Total_IgG1_copy_sort['onecopy'] >= zoom_range_2[0]) & (mutated_Total_IgG1_copy_sort['onecopy'] <= zoom_range_2[1])]
            
            # Plot for DataFrame 1 (Shared IgE) - Zoomed-in around the mean
            sns.scatterplot(data=Mutated_IgE_zoom, x=Mutated_IgE_zoom.index, y='onecopy', color='blue', s=100, ax=axes[1, 0], label='Mutated IgE', edgecolor='w')
            axes[1, 0].axhline(mean1, color='orange', linestyle='--', label=f'IgE Mean: {mean1:.2f}')
            axes[1, 0].axhline(mean2, color='green', linestyle='--', label=f'IgG1 Mean: {mean2:.2f}')
            axes[1, 0].set_title('Zoomed Clonal Size Distribution - Mutated IgE', fontsize=14)
            axes[1, 0].set_xlabel('Clonal Index', fontsize=12)
            axes[1, 0].set_ylabel('Clonal Size', fontsize=12)
            axes[1, 0].legend()
            
            # Plot for DataFrame 2 (Shared IgG1) - Zoomed-in around the mean
            sns.scatterplot(data=Mutated_IgG1_zoom, x=Mutated_IgG1_zoom.index, y='onecopy', color='red', s=100, ax=axes[1, 1], label='Mutated IgG1', edgecolor='w')
            axes[1, 1].axhline(mean1, color='orange', linestyle='--', label=f'IgE Mean: {mean1:.2f}')
            axes[1, 1].axhline(mean2, color='green', linestyle='--', label=f'IgG1 Mean: {mean2:.2f}')
            axes[1, 1].set_title('Zoomed Clonal Size Distribution - Mutated IgG1', fontsize=14)
            axes[1, 1].set_xlabel('Clonal Index', fontsize=12)
            axes[1, 1].set_ylabel('Clonal Size', fontsize=12)
            axes[1, 1].legend()
            
            # Adjust layout to avoid overlap
            plt.tight_layout()
            
            downloads_folder = os.path.join(os.path.expanduser("~"), "Downloads")
            
            plots_folder = os.path.join(downloads_folder, "BCR_Analysis_Plots")
            if not os.path.exists(plots_folder):
                os.makedirs(plots_folder)  # Create the folder if it doesn't exist
                
            plot_path = os.path.join(plots_folder, "clonal size distribution of Mutated IgE and IgG1.png")
            plt.savefig(plot_path, dpi=1200, bbox_inches='tight')
            # Show the plot
            
            
            
            
            #Filter out the Total, mutated and unmutated shared IgE and IgG1
            # Add a column to indicate if 'Sequence.ID' is in the productive shared IgE VDJ and shared IgG1 VDJ
            
            stats_IgE_mut_3['Match_in_shared_IgE_seqid'] = stats_IgE_mut_3['Sequence ID'].isin(set(shared_IgE['seqid']))

            # Check for the correct column name in stats_IgG1_mut_3
            if 'Sequence.ID' in stats_IgG1_mut_3.columns:
                column_seqid = 'Sequence.ID'
            elif 'Sequence ID' in stats_IgG1_mut_3.columns:
                column_seqid = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present in the DataFrame.")
            
            # Perform the analysis by checking if sequence IDs are in the shared_IgG1 set
            stats_IgG1_mut_3['Match_in_shared_IgG1_seqid'] = stats_IgG1_mut_3[column_seqid].isin(set(shared_IgG1['seqid']))
            
            # Filter shared IgE directly using boolean masks
            stats_tab_Shared_IgE = stats_IgE_mut_3[stats_IgE_mut_3['Match_in_shared_IgE_seqid']]
            
            
            # Filter shared IgG1 directly using boolean masks
            stats_tab_Shared_IgG1 = stats_IgG1_mut_3[stats_IgG1_mut_3['Match_in_shared_IgG1_seqid']]
            
            
            # Filter mutated and unmutated shared IgE directly using boolean masks
            stats_tab_mutated_Shared_IgE = stats_IgE_mut_3[
            (stats_IgE_mut_3['Match_in_shared_IgE_seqid']) & (stats_IgE_mut_3['V-REGION Nb of mutations-ex'] > 0)]
            
            stats_tab_unmutated_Shared_IgE = stats_IgE_mut_3[
            (stats_IgE_mut_3['Match_in_shared_IgE_seqid']) & (stats_IgE_mut_3['V-REGION Nb of mutations-ex'] == 0)]
            
            
            # Filter mutated and unmutated shared IgG1 directly using boolean masks
            stats_tab_mutated_Shared_IgG1 = stats_IgG1_mut_3[
            (stats_IgG1_mut_3['Match_in_shared_IgG1_seqid']) & (stats_IgG1_mut_3['V-REGION Nb of mutations-ex'] > 0)]
            
            stats_tab_unmutated_Shared_IgG1 = stats_IgG1_mut_3[
            (stats_IgG1_mut_3['Match_in_shared_IgG1_seqid']) & (stats_IgG1_mut_3['V-REGION Nb of mutations-ex'] == 0)]
            
            
            # Create sets of IgE seqids for faster lookup
            mutated_seqids_Shared_IgE = set(stats_tab_mutated_Shared_IgE['Sequence ID'])
            unmutated_seqids_Shared_IgE = set(stats_tab_unmutated_Shared_IgE['Sequence ID'])
            
    

            # Check for the correct column name in stats_tab_mutated_Shared_IgG1
            if 'Sequence.ID' in stats_tab_mutated_Shared_IgG1.columns:
                column_seqid = 'Sequence.ID'
            elif 'Sequence ID' in stats_tab_mutated_Shared_IgG1.columns:
                column_seqid = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present in the DataFrame.")
            
            # Create a set of unique mutated sequence IDs from the detected column
            mutated_seqids_Shared_IgG1 = set(stats_tab_mutated_Shared_IgG1[column_seqid])


            # Check for the correct column name in stats_tab_unmutated_Shared_IgG1
            if 'Sequence.ID' in stats_tab_unmutated_Shared_IgG1.columns:
                column_seqid = 'Sequence.ID'
            elif 'Sequence ID' in stats_tab_unmutated_Shared_IgG1.columns:
                column_seqid = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present in the DataFrame.")
            
            # Create a set of unique unmutated sequence IDs from the detected column
            unmutated_seqids_Shared_IgG1 = set(stats_tab_unmutated_Shared_IgG1[column_seqid])

            
            # Create a copy of the dataframe
            stats_Shared_IgE_VDJ_seqid_compare = stats_IgE_productive_VDJ.copy()
            stats_Shared_IgG1_VDJ_seqid_compare = stats_IgG1_productive_VDJ.copy()
            
            # Update the comparison IgE DataFrame with boolean masks
            stats_Shared_IgE_VDJ_seqid_compare['Match_in_mutated_shared_IgE_seqid'] = stats_Shared_IgE_VDJ_seqid_compare['seqid'].isin(mutated_seqids_Shared_IgE)
            
            stats_Shared_IgE_VDJ_seqid_compare['Match_in_unmutated_shared_IgE_seqid'] = stats_Shared_IgE_VDJ_seqid_compare['seqid'].isin(unmutated_seqids_Shared_IgE)
            
            
            # Update the comparison IgG1 DataFrame with boolean masks
            stats_Shared_IgG1_VDJ_seqid_compare['Match_in_mutated_shared_IgG1_seqid'] = stats_Shared_IgG1_VDJ_seqid_compare['seqid'].isin(mutated_seqids_Shared_IgG1)
            
            stats_Shared_IgG1_VDJ_seqid_compare['Match_in_unmutated_shared_IgG1_seqid'] = stats_Shared_IgG1_VDJ_seqid_compare['seqid'].isin(unmutated_seqids_Shared_IgG1)
            
            
            
            # Ensure 'Match_in_mutated and unmutated shared IgE seqid' is a string and remove leading/trailing spaces
            stats_Shared_IgE_VDJ_seqid_compare['Match_in_mutated_shared_IgE_seqid'] = stats_Shared_IgE_VDJ_seqid_compare['Match_in_mutated_shared_IgE_seqid'].astype(str).str.strip()
            stats_Shared_IgE_VDJ_seqid_compare['Match_in_unmutated_shared_IgE_seqid'] =stats_Shared_IgE_VDJ_seqid_compare['Match_in_unmutated_shared_IgE_seqid'].astype(str).str.strip()
            
            # Ensure 'Match_in_mutated and unmutated shared IgG1 seqid' is a string and remove leading/trailing spaces
            stats_Shared_IgG1_VDJ_seqid_compare['Match_in_mutated_shared_IgG1_seqid'] = stats_Shared_IgG1_VDJ_seqid_compare['Match_in_mutated_shared_IgG1_seqid'].astype(str).str.strip()
            stats_Shared_IgG1_VDJ_seqid_compare['Match_in_unmutated_shared_IgG1_seqid'] = stats_Shared_IgG1_VDJ_seqid_compare['Match_in_unmutated_shared_IgG1_seqid'].astype(str).str.strip()
            
            # Now filter the mutated and unmutated shared IgE clones
            stats_tab_mutated_shared_IgE_VDJ = stats_Shared_IgE_VDJ_seqid_compare[(stats_Shared_IgE_VDJ_seqid_compare['Match_in_mutated_shared_IgE_seqid'] == 'True')]
            stats_tab_unmutated_shared_IgE_VDJ = stats_Shared_IgE_VDJ_seqid_compare[(stats_Shared_IgE_VDJ_seqid_compare['Match_in_unmutated_shared_IgE_seqid'] == 'True')]
            
            # Now filter the mutated and unmutated shared IgG clones
            stats_tab_mutated_shared_IgG1_VDJ = stats_Shared_IgG1_VDJ_seqid_compare[(stats_Shared_IgG1_VDJ_seqid_compare['Match_in_mutated_shared_IgG1_seqid'] == 'True')]
            stats_tab_unmutated_shared_IgG1_VDJ = stats_Shared_IgG1_VDJ_seqid_compare[(stats_Shared_IgG1_VDJ_seqid_compare['Match_in_unmutated_shared_IgG1_seqid'] == 'True')]
            
            
            
            # Plot copy number distribution of mutated and Unmutated Shared IgE and IgG1
            # Sort the counts in ascending order
            mutated_shared_IgE_copy_sort = stats_tab_mutated_shared_IgE_VDJ.sort_values(by='onecopy', ascending=True).reset_index(drop=True)
            mutated_shared_IgG1_copy_sort = stats_tab_mutated_shared_IgG1_VDJ.sort_values(by='onecopy', ascending=True).reset_index(drop=True)
            
            # Calculate the mean and median values for both dataframes
            mean1 = mutated_shared_IgE_copy_sort ['onecopy'].mean()
            mean2 = mutated_shared_IgG1_copy_sort ['onecopy'].mean()
            
            # Set Seaborn style
            sns.set(style="whitegrid")
            
            # Create the plot with two rows and two columns of subplots (original plots and zoomed-in plots)
            fig8, axes = plt.subplots(2, 2, figsize=(14, 12))
            
            # Plot for DataFrame 1 (mutated Shared IgE) - Original
            sns.scatterplot(data=mutated_shared_IgE_copy_sort, x=mutated_shared_IgE_copy_sort.index, y='onecopy', color='blue', s=100, ax=axes[0, 0], label='Mutated shared IgE', edgecolor='w')
            axes[0, 0].axhline(mean1, color='orange', linestyle='--', label=f'IgE Mean: {mean1:.2f}')
            axes[0, 0].axhline(mean2, color='green', linestyle='--', label=f'IgG1 Mean: {mean2:.2f}')
            axes[0, 0].set_title('Clonal Size Distribution - Mutated shared IgE', fontsize=14)
            axes[0, 0].set_xlabel('Clonal Index', fontsize=12)
            axes[0, 0].set_ylabel('Clonal Size', fontsize=12)
            axes[0, 0].legend()
            
            # Plot for DataFrame 2 (mutated Shared IgG1) - Original
            sns.scatterplot(data=mutated_shared_IgG1_copy_sort, x=mutated_shared_IgG1_copy_sort.index, y='onecopy', color='red', s=100, ax=axes[0, 1], label='Mutated shared IgG1', edgecolor='w')
            axes[0, 1].axhline(mean1, color='orange', linestyle='--', label=f'IgE Mean: {mean1:.2f}')
            axes[0, 1].axhline(mean2, color='green', linestyle='--', label=f'IgG1 Mean: {mean2:.2f}')
            axes[0, 1].set_title('Clonal Size Distribution - Mutated shared IgG1', fontsize=14)
            axes[0, 1].set_xlabel('Clonal Index', fontsize=12)
            axes[0, 1].set_ylabel('Clonal Size', fontsize=12)
            axes[0, 1].legend()
            
            # Define a range around the mean for zooming in (e.g., +-10 around the mean)
            zoom_range_1 = (mean2 - 20, mean2 + 20)
            zoom_range_2 = (mean2 - 20, mean2 + 20)
            
            # Filter DataFrames to show values around the mean
            Mutated_IgE_zoom = mutated_shared_IgE_copy_sort[(mutated_shared_IgE_copy_sort['onecopy'] >= zoom_range_1[0]) & (mutated_shared_IgE_copy_sort['onecopy'] <= zoom_range_1[1])]
            Mutated_IgG1_zoom = mutated_shared_IgG1_copy_sort[(mutated_shared_IgG1_copy_sort['onecopy'] >= zoom_range_2[0]) & (mutated_shared_IgG1_copy_sort['onecopy'] <= zoom_range_2[1])]
            
            # Plot for DataFrame 1 (mutated Shared IgE) - Zoomed-in around the mean
            sns.scatterplot(data=Mutated_IgE_zoom, x=Mutated_IgE_zoom.index, y='onecopy', color='blue', s=100, ax=axes[1, 0], label='Mutated shared IgE', edgecolor='w')
            axes[1, 0].axhline(mean1, color='orange', linestyle='--', label=f'IgE Mean: {mean1:.2f}')
            axes[1, 0].axhline(mean2, color='green', linestyle='--', label=f'IgG1 Mean: {mean2:.2f}')
            axes[1, 0].set_title('Zoomed Clonal Size Distribution - Mutated shared IgE', fontsize=14)
            axes[1, 0].set_xlabel('Clonal Index', fontsize=12)
            axes[1, 0].set_ylabel('Clonal Size', fontsize=12)
            axes[1, 0].legend()
            
            # Plot for DataFrame 2 (mutated Shared IgG1) - Zoomed-in around the mean
            sns.scatterplot(data=Mutated_IgG1_zoom, x=Mutated_IgG1_zoom.index, y='onecopy', color='red', s=100, ax=axes[1, 1], label='Mutated shared IgG1', edgecolor='w')
            axes[1, 1].axhline(mean1, color='orange', linestyle='--', label=f'IgE Mean: {mean1:.2f}')
            axes[1, 1].axhline(mean2, color='green', linestyle='--', label=f'IgG1 Mean: {mean2:.2f}')
            axes[1, 1].set_title('Zoomed Clonal Size Distribution - Mutated shared IgG1', fontsize=14)
            axes[1, 1].set_xlabel('Clonal Index', fontsize=12)
            axes[1, 1].set_ylabel('Clonal Size', fontsize=12)
            axes[1, 1].legend()
            
            # Adjust layout to avoid overlap
            plt.tight_layout()
            
            downloads_folder = os.path.join(os.path.expanduser("~"), "Downloads")
            
            plots_folder = os.path.join(downloads_folder, "BCR_Analysis_Plots")
            if not os.path.exists(plots_folder):
                os.makedirs(plots_folder)  # Create the folder if it doesn't exist
                
            plot_path = os.path.join(plots_folder, "Clonal size distribution of mutated shared IgE and IgG1.png")
            plt.savefig(plot_path, dpi=1200, bbox_inches='tight')
            
            # Show the plot
            
            
            #Plot graph to show the mutation distribution of Shared IgE and IgG1
            from collections import Counter
            
            # Data: Two sets of mutation counts
            mutation_counts_IgE = [stats_tab_Shared_IgE['V-REGION Nb of mutations-ex']]
            
            mutation_counts_IgG1 = [stats_tab_Shared_IgG1['V-REGION Nb of mutations-ex']]
            
            #Convert the Series into a list by using .tolist()
            mutation_counts_IgE = stats_tab_Shared_IgE['V-REGION Nb of mutations-ex'].tolist()
            mutation_counts_IgG1 = stats_tab_Shared_IgG1['V-REGION Nb of mutations-ex'].tolist()
            
            # Count the occurrences of each mutation in both datasets
            mutation_distribution_IgE = Counter(mutation_counts_IgE)
            mutation_distribution_IgG1 = Counter(mutation_counts_IgG1)
            
            # Sort the data by mutation number in ascending order for both distributions
            sorted_mutations_IgE = sorted(mutation_distribution_IgE.items())
            sorted_mutations_IgG1 = sorted(mutation_distribution_IgG1.items())
            
            # Separate into X and Y values for both distributions
            x_values_IgE = [mutation for mutation, count in sorted_mutations_IgE]
            y_values_IgE = [count for mutation, count in sorted_mutations_IgE]
            
            x_values_IgG1 = [mutation for mutation, count in sorted_mutations_IgG1]
            y_values_IgG1 = [count for mutation, count in sorted_mutations_IgG1]
            
            # Check if lengths match before plotting
            if len(x_values_IgE) != len(y_values_IgE):
                raise ValueError(f"Mismatch in lengths of x_values_IgE ({len(x_values_IgE)}) and y_values_IgE ({len(y_values_IgE)})")
            if len(x_values_IgG1) != len(y_values_IgG1):
                raise ValueError(f"Mismatch in lengths of x_values_IgG1 ({len(x_values_IgG1)}) and y_values_IgG1 ({len(y_values_IgG1)})")
            
            # Plot both mutation distributions as dot plots with different colors and markers
            fig9 = plt.figure(figsize=(10, 6))
            
            # Plot the first distribution in blue with circle markers and connect the dots with a thin line
            plt.plot(x_values_IgE, y_values_IgE, 'bo-', markersize=8, label='Shared IgE', linewidth=1)
            
            # Plot the second distribution in red with square markers and connect the dots with a thin line
            plt.plot(x_values_IgG1, y_values_IgG1, 'rs-', markersize=8, label='Shared IgG1', linewidth=1)
            
            # Add vertical grid lines for mutation categories
            plt.axvline(x=5, color='green', linestyle='--', linewidth=1)
            plt.axvline(x=11, color='orange', linestyle='--', linewidth=1)
            plt.axvline(x=12, color='red', linestyle='--', linewidth=1)
            
            # Add labels for mutation categories in their respective ranges
            plt.text(2.5, max(y_values_IgE + y_values_IgG1), 'Low Mutation (0-5)', color='green', fontsize=10.5, ha='right')
            plt.text(8.5, max(y_values_IgE + y_values_IgG1), 'Moderate Mutation (6-11)', color='orange', fontsize=10.5, ha='center')
            plt.text(13, max(y_values_IgE + y_values_IgG1), 'High Mutation (>=12)', color='red', fontsize=10.5, ha='left')
            
            # Add labels, title, and legend
            plt.xlabel('VH region number of mutations', fontsize=12)
            plt.ylabel('Count', fontsize=12)
            plt.title('Mutation Distributions', fontsize=14)
            plt.legend()
            
            # Show grid for better readability
            plt.grid(axis='both', linestyle='--', alpha=0.2)
            
            # Display the plot
            plt.tight_layout()
            downloads_folder = os.path.join(os.path.expanduser("~"), "Downloads")
            
            plots_folder = os.path.join(downloads_folder, "BCR_Analysis_Plots")
            if not os.path.exists(plots_folder):
                os.makedirs(plots_folder)  # Create the folder if it doesn't exist
                
            plot_path = os.path.join(plots_folder, "Mutation distribution of shared IgE and IgG1.png")
            plt.savefig(plot_path, dpi=1200, bbox_inches='tight')
            
            
            
            
            # Filter out the low, moderate and highly mutated shared IgE and IgG1 mutation table
            # Filter out the low mutated shared IgE
            stats_tab_mutated_low_mut_Shared_IgE = stats_IgE_mut_3[
            (stats_IgE_mut_3['Match_in_shared_IgE_seqid']) &
            (stats_IgE_mut_3['V-REGION Nb of mutations-ex'] >= 0) &
            (stats_IgE_mut_3['V-REGION Nb of mutations-ex'] <= 5)]
            
            # Filter out the low mutated shared IgG1
            stats_tab_mutated_low_mut_Shared_IgG1 = stats_IgG1_mut_3[
            (stats_IgG1_mut_3['Match_in_shared_IgG1_seqid']) &
            (stats_IgG1_mut_3['V-REGION Nb of mutations-ex'] >= 0) &
            (stats_IgG1_mut_3['V-REGION Nb of mutations-ex'] <= 5)]
            
            # Filter out the moderate mutated shared IgE
            stats_tab_mutated_mod_mut_Shared_IgE = stats_IgE_mut_3[
            (stats_IgE_mut_3['Match_in_shared_IgE_seqid']) &
            (stats_IgE_mut_3['V-REGION Nb of mutations-ex'] >= 6) &
            (stats_IgE_mut_3['V-REGION Nb of mutations-ex'] <= 11)]
            
            # Filter out the moderate mutated shared IgG1
            stats_tab_mutated_mod_mut_Shared_IgG1 = stats_IgG1_mut_3[
            (stats_IgG1_mut_3['Match_in_shared_IgG1_seqid']) &
            (stats_IgG1_mut_3['V-REGION Nb of mutations-ex'] >= 6) &
            (stats_IgG1_mut_3['V-REGION Nb of mutations-ex'] <= 11)]
            
            # Filter out the highly mutated shared IgE
            stats_tab_mutated_high_mut_Shared_IgE = stats_IgE_mut_3[
            (stats_IgE_mut_3['Match_in_shared_IgE_seqid']) &
            (stats_IgE_mut_3['V-REGION Nb of mutations-ex'] >= 12)]
            
            # Filter out the highly mutated shared IgG1
            stats_tab_mutated_high_mut_Shared_IgG1 = stats_IgG1_mut_3[
            (stats_IgG1_mut_3['Match_in_shared_IgG1_seqid']) &
            (stats_IgG1_mut_3['V-REGION Nb of mutations-ex'] >= 12)]
            
            
            
            # Create sets of highly, moderate and low mutated IgE seqids for faster lookup
            low_mutated_seqids_Shared_IgE = set(stats_tab_mutated_low_mut_Shared_IgE['Sequence ID'])
            moderate_mutated_seqids_Shared_IgE = set(stats_tab_mutated_mod_mut_Shared_IgE['Sequence ID'])
            High_mutated_seqids_Shared_IgE = set(stats_tab_mutated_high_mut_Shared_IgE['Sequence ID'])
            
            
            # Create sets of highly, moderate and low mutated IgG1 seqids for faster lookup
            # Check for the correct column name in stats_tab_mutated_low_mut_Shared_IgG1
            if 'Sequence.ID' in stats_tab_mutated_low_mut_Shared_IgG1.columns:
                column_seqid = 'Sequence.ID'
            elif 'Sequence ID' in stats_tab_mutated_low_mut_Shared_IgG1.columns:
                column_seqid = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present in the DataFrame.")
            
            # Create a set of unique low mutated sequence IDs from the detected column
            low_mutated_seqids_Shared_IgG1 = set(stats_tab_mutated_low_mut_Shared_IgG1[column_seqid])



            # Check for the correct column name in stats_tab_mutated_mod_mut_Shared_IgG1
            if 'Sequence.ID' in stats_tab_mutated_mod_mut_Shared_IgG1.columns:
                column_seqid = 'Sequence.ID'
            elif 'Sequence ID' in stats_tab_mutated_mod_mut_Shared_IgG1.columns:
                column_seqid = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present in the DataFrame.")
            
            # Create a set of unique moderate mutated sequence IDs from the detected column
            moderate_mutated_seqids_Shared_IgG1 = set(stats_tab_mutated_mod_mut_Shared_IgG1[column_seqid])


            # Check for the correct column name in stats_tab_mutated_high_mut_Shared_IgG1
            if 'Sequence.ID' in stats_tab_mutated_high_mut_Shared_IgG1.columns:
                column_seqid = 'Sequence.ID'
            elif 'Sequence ID' in stats_tab_mutated_high_mut_Shared_IgG1.columns:
                column_seqid = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present in the DataFrame.")
            
            # Create a set of unique high mutated sequence IDs from the detected column
            High_mutated_seqids_Shared_IgG1 = set(stats_tab_mutated_high_mut_Shared_IgG1[column_seqid])

            
            # Create a copy of the dataframe
            stats_mut_Shared_IgE_VDJ_seqid_compare = stats_IgE_productive_VDJ.copy()
            stats_mut_Shared_IgG1_VDJ_seqid_compare = stats_IgG1_productive_VDJ.copy()
            
            # Update the comparison IgE DataFrame with boolean masks
            stats_mut_Shared_IgE_VDJ_seqid_compare['Match_in_low_mutated_shared_IgE_seqid'] = stats_mut_Shared_IgE_VDJ_seqid_compare['seqid'].isin(low_mutated_seqids_Shared_IgE)
            stats_mut_Shared_IgE_VDJ_seqid_compare['Match_in_moderate_mutated_shared_IgE_seqid'] = stats_mut_Shared_IgE_VDJ_seqid_compare['seqid'].isin(moderate_mutated_seqids_Shared_IgE)
            stats_mut_Shared_IgE_VDJ_seqid_compare['Match_in_high_mutated_shared_IgE_seqid'] = stats_mut_Shared_IgE_VDJ_seqid_compare['seqid'].isin(High_mutated_seqids_Shared_IgE)
            
            
            # Update the comparison IgG1 DataFrame with boolean masks
            stats_mut_Shared_IgG1_VDJ_seqid_compare['Match_in_low_mutated_shared_IgG1_seqid'] = stats_mut_Shared_IgG1_VDJ_seqid_compare['seqid'].isin(low_mutated_seqids_Shared_IgG1)
            stats_mut_Shared_IgG1_VDJ_seqid_compare['Match_in_moderate_mutated_shared_IgG1_seqid'] = stats_mut_Shared_IgG1_VDJ_seqid_compare['seqid'].isin(moderate_mutated_seqids_Shared_IgG1)
            stats_mut_Shared_IgG1_VDJ_seqid_compare['Match_in_high_mutated_shared_IgG1_seqid'] = stats_mut_Shared_IgG1_VDJ_seqid_compare['seqid'].isin(High_mutated_seqids_Shared_IgG1)
            
            
            # Ensure 'Match_in_low_moderate and high mutated_shared_IgE_seqid' is a string and remove leading/trailing spaces
            stats_mut_Shared_IgE_VDJ_seqid_compare['Match_in_low_mutated_shared_IgE_seqid'] = stats_mut_Shared_IgE_VDJ_seqid_compare['Match_in_low_mutated_shared_IgE_seqid'].astype(str).str.strip()
            stats_mut_Shared_IgE_VDJ_seqid_compare['Match_in_moderate_mutated_shared_IgE_seqid'] =stats_mut_Shared_IgE_VDJ_seqid_compare['Match_in_moderate_mutated_shared_IgE_seqid'].astype(str).str.strip()
            stats_mut_Shared_IgE_VDJ_seqid_compare['Match_in_high_mutated_shared_IgE_seqid'] = stats_mut_Shared_IgE_VDJ_seqid_compare['Match_in_high_mutated_shared_IgE_seqid'].astype(str).str.strip()
            
            # Ensure 'Match_in_low_moderate and high mutated_shared_IgG1_seqid' is a string and remove leading/trailing spaces
            stats_mut_Shared_IgG1_VDJ_seqid_compare['Match_in_low_mutated_shared_IgG1_seqid'] = stats_mut_Shared_IgG1_VDJ_seqid_compare['Match_in_low_mutated_shared_IgG1_seqid'].astype(str).str.strip()
            stats_mut_Shared_IgG1_VDJ_seqid_compare['Match_in_moderate_mutated_shared_IgG1_seqid'] = stats_mut_Shared_IgG1_VDJ_seqid_compare['Match_in_moderate_mutated_shared_IgG1_seqid'].astype(str).str.strip()
            stats_mut_Shared_IgG1_VDJ_seqid_compare['Match_in_high_mutated_shared_IgG1_seqid'] = stats_mut_Shared_IgG1_VDJ_seqid_compare['Match_in_high_mutated_shared_IgG1_seqid'].astype(str).str.strip()
            
            
            # Now filter the low_moderate and high mutated_shared_IgE clones
            stats_tab_low_mutated_shared_IgE_VDJ = stats_mut_Shared_IgE_VDJ_seqid_compare[(stats_mut_Shared_IgE_VDJ_seqid_compare['Match_in_low_mutated_shared_IgE_seqid'] == 'True')]
            stats_tab_moderate_mutated_shared_IgE_VDJ = stats_mut_Shared_IgE_VDJ_seqid_compare[(stats_mut_Shared_IgE_VDJ_seqid_compare['Match_in_moderate_mutated_shared_IgE_seqid'] == 'True')]
            stats_tab_high_mutated_shared_IgE_VDJ = stats_mut_Shared_IgE_VDJ_seqid_compare[(stats_mut_Shared_IgE_VDJ_seqid_compare['Match_in_high_mutated_shared_IgE_seqid'] == 'True')]
            
            
            ## Now filter the low_moderate and high mutated_shared_IgG1 clones
            stats_tab_low_mutated_shared_IgG1_VDJ = stats_mut_Shared_IgG1_VDJ_seqid_compare[(stats_mut_Shared_IgG1_VDJ_seqid_compare['Match_in_low_mutated_shared_IgG1_seqid'] == 'True')]
            stats_tab_moderate_mutated_shared_IgG1_VDJ = stats_mut_Shared_IgG1_VDJ_seqid_compare[(stats_mut_Shared_IgG1_VDJ_seqid_compare['Match_in_moderate_mutated_shared_IgG1_seqid'] == 'True')]
            stats_tab_high_mutated_shared_IgG1_VDJ = stats_mut_Shared_IgG1_VDJ_seqid_compare[(stats_mut_Shared_IgG1_VDJ_seqid_compare['Match_in_high_mutated_shared_IgG1_seqid'] == 'True')]
            
            
            
            
            
            
            # plot the low, moderate and highly mutated shared IgE and IgG1 clonal size distribution
            # Sort the counts in ascending order
            stats_tab_low_mutated_shared_IgE_VDJ_sort = stats_tab_low_mutated_shared_IgE_VDJ.sort_values(by='onecopy', ascending=True).reset_index(drop=True)
            stats_tab_low_mutated_shared_IgG1_VDJ_sort = stats_tab_low_mutated_shared_IgG1_VDJ.sort_values(by='onecopy', ascending=True).reset_index(drop=True)
            
            stats_tab_moderate_mutated_shared_IgE_VDJ_sort = stats_tab_moderate_mutated_shared_IgE_VDJ.sort_values(by='onecopy', ascending=True).reset_index(drop=True)
            stats_tab_moderate_mutated_shared_IgG1_VDJ_sort = stats_tab_moderate_mutated_shared_IgG1_VDJ.sort_values(by='onecopy', ascending=True).reset_index(drop=True)
            
            stats_tab_high_mutated_shared_IgE_VDJ_sort = stats_tab_high_mutated_shared_IgE_VDJ.sort_values(by='onecopy', ascending=True).reset_index(drop=True)
            stats_tab_high_mutated_shared_IgG1_VDJ_sort = stats_tab_high_mutated_shared_IgG1_VDJ.sort_values(by='onecopy', ascending=True).reset_index(drop=True)
            
            
            # Calculate the mean and median values for both dataframes
            mean1 = stats_tab_low_mutated_shared_IgE_VDJ_sort ['onecopy'].mean()
            mean2 = stats_tab_low_mutated_shared_IgG1_VDJ_sort ['onecopy'].mean()
            
            mean3 = stats_tab_moderate_mutated_shared_IgE_VDJ_sort ['onecopy'].mean()
            mean4 = stats_tab_moderate_mutated_shared_IgG1_VDJ_sort ['onecopy'].mean()
            
            mean5 = stats_tab_high_mutated_shared_IgE_VDJ_sort ['onecopy'].mean()
            mean6 = stats_tab_high_mutated_shared_IgG1_VDJ_sort ['onecopy'].mean()
            
            
            # Set Seaborn style
            sns.set(style="whitegrid")
            
            # Create the plot with two rows and two columns of subplots (original plots and zoomed-in plots)
            fig10, axes = plt.subplots(3, 2, figsize=(14, 12))
            
            # Plot for DataFrame 1 (low mutated Shared IgE) - Original
            sns.scatterplot(data=stats_tab_low_mutated_shared_IgE_VDJ_sort, x=stats_tab_low_mutated_shared_IgE_VDJ_sort.index, y='onecopy', color='blue', s=100, ax=axes[0, 0], label='low mutated shared IgE', edgecolor='w')
            axes[0, 0].axhline(mean1, color='orange', linestyle='--', label=f'mean low mutated IgE: {mean1:.2f}')
            axes[0, 0].axhline(mean2, color='green', linestyle='--', label=f'mean low mutated IgG1: {mean2:.2f}')
            axes[0, 0].set_title('low mutated shared IgE', fontsize=14)
            axes[0, 0].set_xlabel('Clonal Index', fontsize=12)
            axes[0, 0].set_ylabel('Clonal Size', fontsize=12)
            axes[0, 0].legend()
            
            # Plot for DataFrame 2 (low mutated Shared IgG1) - Original
            sns.scatterplot(data=stats_tab_low_mutated_shared_IgG1_VDJ_sort, x=stats_tab_low_mutated_shared_IgG1_VDJ_sort.index, y='onecopy', color='red', s=100, ax=axes[0, 1], label='low mutated shared IgG1', edgecolor='w')
            axes[0, 1].axhline(mean1, color='orange', linestyle='--', label=f'mean low mutated IgE: {mean1:.2f}')
            axes[0, 1].axhline(mean2, color='green', linestyle='--', label=f'mean low mutated IgG1: {mean2:.2f}')
            axes[0, 1].set_title('low mutated shared IgG1', fontsize=14)
            axes[0, 1].set_xlabel('Clonal Index', fontsize=12)
            axes[0, 1].set_ylabel('Clonal Size', fontsize=12)
            axes[0, 1].legend()
            
            # Plot for DataFrame 3 (moderate mutated Shared IgE) -
            sns.scatterplot(data=stats_tab_moderate_mutated_shared_IgE_VDJ_sort, x=stats_tab_moderate_mutated_shared_IgE_VDJ_sort.index, y='onecopy', color='blue', s=100, ax=axes[1, 0], label='moderate mutated shared IgE', edgecolor='w')
            axes[1, 0].axhline(mean3, color='orange', linestyle='--', label=f'mean moderate mutated IgE: {mean3:.2f}')
            axes[1, 0].axhline(mean4, color='green', linestyle='--', label=f'mean moderate mutated IgG1: {mean4:.2f}')
            axes[1, 0].set_title('moderate mutated shared IgE', fontsize=14)
            axes[1, 0].set_xlabel('Clonal Index', fontsize=12)
            axes[1, 0].set_ylabel('Clonal Size', fontsize=12)
            axes[1, 0].legend()
            
            # Plot for DataFrame 4 (moderate mutated Shared IgG1) -
            sns.scatterplot(data=stats_tab_moderate_mutated_shared_IgG1_VDJ_sort, x=stats_tab_moderate_mutated_shared_IgG1_VDJ_sort.index, y='onecopy', color='red', s=100, ax=axes[1, 1], label='moderate mutated shared IgG1', edgecolor='w')
            axes[1, 1].axhline(mean3, color='orange', linestyle='--', label=f'mean moderate mutated IgE: {mean3:.2f}')
            axes[1, 1].axhline(mean4, color='green', linestyle='--', label=f'mean moderate mutated IgG1: {mean4:.2f}')
            axes[1, 1].set_title('moderate mutated shared IgG1', fontsize=14)
            axes[1, 1].set_xlabel('Clonal Index', fontsize=12)
            axes[1, 1].set_ylabel('Clonal Size', fontsize=12)
            axes[1, 1].legend()
            
            
            # Plot for DataFrame 5 (high mutated Shared IgE) -
            sns.scatterplot(data=stats_tab_high_mutated_shared_IgE_VDJ_sort, x=stats_tab_high_mutated_shared_IgE_VDJ_sort.index, y='onecopy', color='blue', s=100, ax=axes[2, 0], label='high mutated shared IgE', edgecolor='w')
            axes[2, 0].axhline(mean5, color='orange', linestyle='--', label=f'mean high mutated IgE: {mean5:.2f}')
            axes[2, 0].axhline(mean6, color='green', linestyle='--', label=f'mean high mutated IgG1: {mean6:.2f}')
            axes[2, 0].set_title('high mutated shared IgE', fontsize=14)
            axes[2, 0].set_xlabel('Clonal Index', fontsize=12)
            axes[2, 0].set_ylabel('Clonal Size', fontsize=12)
            axes[2, 0].legend()
            
            
            # Plot for DataFrame 6 (high mutated Shared IgG1) -
            sns.scatterplot(data=stats_tab_high_mutated_shared_IgG1_VDJ_sort, x=stats_tab_high_mutated_shared_IgG1_VDJ_sort.index, y='onecopy', color='red', s=100, ax=axes[2, 1], label='high mutated shared IgG1', edgecolor='w')
            axes[2, 1].axhline(mean5, color='orange', linestyle='--', label=f'mean high mutated IgE: {mean5:.2f}')
            axes[2, 1].axhline(mean6, color='green', linestyle='--', label=f'mean high mutated IgG1: {mean6:.2f}')
            axes[2, 1].set_title('high mutated shared IgG1', fontsize=14)
            axes[2, 1].set_xlabel('Clonal Index', fontsize=12)
            axes[2, 1].set_ylabel('Clonal Size', fontsize=12)
            axes[2, 1].legend()
            
            
            # Adjust layout to avoid overlap
            plt.tight_layout()
            downloads_folder = os.path.join(os.path.expanduser("~"), "Downloads")
            
            plots_folder = os.path.join(downloads_folder, "BCR_Analysis_Plots")
            if not os.path.exists(plots_folder):
                os.makedirs(plots_folder)  # Create the folder if it doesn't exist
                
            plot_path = os.path.join(plots_folder, "Clonal size distribution of low, moderate and high mutated shared IgE and IgG1.png")
            plt.savefig(plot_path, dpi=1200, bbox_inches='tight')
            
            # Show the plot
            
            
            
            
            # Determine the mean of the low, moderate and high mutated IgE and IG1
            
            print(f'mean copynumber low mutated IgE: {mean1}')
            
            print(f'mean copynumber low mutated IgG1: {mean2}')
            
            print(f'mean copynumber moderate mutated IgE: {mean3}')
            
            print(f'mean copynumber moderate mutated IgG1: {mean4}')
            
            print(f'mean copynumber high mutated IgE: {mean5}')
            
            print(f'mean copynumber high mutated IgG1: {mean6}')
            
            
            # Filter the 10 most expanded low, moderate and highly mutated shared IgE and IgG1
            
            # Define a single reusable filter function
            def filter_top_ten(df, column):
                """
                Sort the DataFrame by the specified column in descending order 
                and return the top 10 rows.
                
                Parameters:
                df (pd.DataFrame): Input DataFrame to filter.
                column (str): Column name to sort by.
                
                Returns:
                pd.DataFrame: Sorted DataFrame with top 10 rows.
                """
                if df is not None and column in df.columns:
                    return df.sort_values(by=column, ascending=False).head(10)
                else:
                    print("Error: Invalid DataFrame or column.")
                    return pd.DataFrame()  # Return an empty DataFrame if input is invalid
            
            # Create a sample DataFrame for testing (replace with actual data)
            # Assuming 'onecopy' is the column to sort by
            sequence_column = 'cdr3aa'  # Replace with actual column name
            
            # Filter low, moderate, and high mutated IgE and IgG1 DataFrames
            stats_tab_low_mutated_shared_IgE_VDJ_sorted_topten = filter_top_ten(stats_tab_low_mutated_shared_IgE_VDJ, 'onecopy')
            stats_tab_moderate_mutated_shared_IgE_VDJ_sorted_topten = filter_top_ten(stats_tab_moderate_mutated_shared_IgE_VDJ, 'onecopy')
            stats_tab_high_mutated_shared_IgE_VDJ_sorted_topten = filter_top_ten(stats_tab_high_mutated_shared_IgE_VDJ, 'onecopy')
            
            stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten = filter_top_ten(stats_tab_low_mutated_shared_IgG1_VDJ, 'onecopy')
            stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten = filter_top_ten(stats_tab_moderate_mutated_shared_IgG1_VDJ, 'onecopy')
            stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten = filter_top_ten(stats_tab_high_mutated_shared_IgG1_VDJ, 'onecopy')
            
            
            # Define a single reusable function for calculating net charge
            def calculate_net_charge(df, sequence_column, pH=7.0, pKscale="Lehninger"):
                """
                Calculate the net charge of peptide sequences in a DataFrame column.
                
                Parameters:
                df (pd.DataFrame): The input DataFrame containing peptide sequences.
                sequence_column (str): The column name with peptide sequences.
                pH (float): The pH at which to compute the net charge. Default is 7.0.
                pKscale (str): The pKa scale to be used for computing the net charge. Default is 'Lehninger'.
                
                Returns:
                pd.Series: A pandas Series with the net charge for each peptide sequence.
                """
                def net_charge(sequence):
                    peptide = Peptide(sequence)
                    return peptide.charge(pH=pH, pKscale=pKscale)
                
                return df[sequence_column].apply(net_charge)
            
            # Assuming filtered DataFrames are already defined
            sequence_column = 'cdr3aa'  # Replace with your actual column name
            
            # Apply net charge calculation to all filtered DataFrames
            stats_tab_low_mutated_shared_IgE_VDJ_sorted_topten['Net charge CDR3aa'] = calculate_net_charge(
                stats_tab_low_mutated_shared_IgE_VDJ_sorted_topten, sequence_column)
            
            stats_tab_moderate_mutated_shared_IgE_VDJ_sorted_topten['Net charge CDR3aa'] = calculate_net_charge(
                stats_tab_moderate_mutated_shared_IgE_VDJ_sorted_topten, sequence_column)
            
            stats_tab_high_mutated_shared_IgE_VDJ_sorted_topten['Net charge CDR3aa'] = calculate_net_charge(
                stats_tab_high_mutated_shared_IgE_VDJ_sorted_topten, sequence_column)
            
            stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten['Net charge CDR3aa'] = calculate_net_charge(
                stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten, sequence_column)
            
            stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten['Net charge CDR3aa'] = calculate_net_charge(
                stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten, sequence_column)
            
            stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten['Net charge CDR3aa'] = calculate_net_charge(
                stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten, sequence_column)
            
            
            # Determine the mean Netcharge of the low, moderate and high mutated IgE
            if stats_tab_low_mutated_shared_IgE_VDJ_sorted_topten is not None:
                mean_stats_tab_low_mutated_shared_IgE_VDJ_sorted_topten = stats_tab_low_mutated_shared_IgE_VDJ_sorted_topten['Net charge CDR3aa'].mean()
                print(f'Mean CDR3 netcharge of low mutated IgE: {mean_stats_tab_low_mutated_shared_IgE_VDJ_sorted_topten}')
            else:
                print("The DataFrame 'stats_tab_low_mutated_shared_IgE_VDJ_sorted_topten' is None")
            
            mean_stats_tab_moderate_mutated_shared_IgE_VDJ_sorted_topten = stats_tab_moderate_mutated_shared_IgE_VDJ_sorted_topten['Net charge CDR3aa'].mean()
            print(f'Mean CDR3 netcharge of moderate mutated IgE: {mean_stats_tab_moderate_mutated_shared_IgE_VDJ_sorted_topten}')
            
            mean_stats_tab_high_mutated_shared_IgE_VDJ_sorted_topten = stats_tab_high_mutated_shared_IgE_VDJ_sorted_topten['Net charge CDR3aa'].mean()
            print(f'Mean CDR3 netcharge of high mutated IgE: {mean_stats_tab_high_mutated_shared_IgE_VDJ_sorted_topten}')
            
            
            # Determine the mean Netcharge of the low, moderate and high mutated IgG1
            
            mean_stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten = stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten['Net charge CDR3aa'].mean()
            print(f'Mean CDR3 netcharge of low mutated IgG1: {mean_stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten}')
            
            mean_stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten = stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten['Net charge CDR3aa'].mean()
            print(f'Mean CDR3 netcharge of moderate mutated IgG1: {mean_stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten}')
            
            mean_stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten = stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten['Net charge CDR3aa'].mean()
            print(f'Mean CDR3 netcharge of high mutated IgG1: {mean_stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten}')
            
            
            
            # Determine the Netcharge ratio of the low, moderate and high mutated IgG1/IgE
            ratio_stats_tab_low_mutated_shared_IgG1_IgE_VDJ_sorted_topten = (mean_stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten /mean_stats_tab_low_mutated_shared_IgE_VDJ_sorted_topten)
            print(f'low mutated IgG1/IgE netcharge ratio: {ratio_stats_tab_low_mutated_shared_IgG1_IgE_VDJ_sorted_topten}')
            
            ratio_stats_tab_moderate_mutated_shared_IgG1_IgE_VDJ_sorted_topten = (mean_stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten /mean_stats_tab_moderate_mutated_shared_IgE_VDJ_sorted_topten)
            print(f'moderate mutated IgG1/IgE netcharge ratio: {ratio_stats_tab_moderate_mutated_shared_IgG1_IgE_VDJ_sorted_topten}')
            
            ratio_stats_tab_high_mutated_shared_IgG1_IgE_VDJ_sorted_topten = (mean_stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten /mean_stats_tab_high_mutated_shared_IgE_VDJ_sorted_topten)
            print(f'high mutated IgG1/IgE netcharge ratio: {ratio_stats_tab_high_mutated_shared_IgG1_IgE_VDJ_sorted_topten}')
            
            
            #Compare the low, moderate and high mutated IgE CDR3aa to the IgG1 Aminoacid tab
            #Filter productive sequences from the IgE Aminoacid table
            stats_IgE_Aminotab_productive = stats_IgE_Aminotab[(stats_IgE_Aminotab['V-DOMAIN Functionality'] == "productive")]
            

            # Check if either 'V.DOMAIN Functionality' or 'V-DOMAIN Functionality' exists in the DataFrame

            if 'V.DOMAIN.Functionality' in stats_IgG1_Aminotab.columns:
                column_name = 'V.DOMAIN.Functionality'
            elif 'V-DOMAIN Functionality' in stats_IgG1_Aminotab.columns:
                column_name = 'V-DOMAIN Functionality'
            else:
                raise KeyError("Neither 'V.DOMAIN.Functionality' nor 'V-DOMAIN Functionality' column is present in the DataFrame.")
            
            # Filter for productive sequences
            stats_IgG1_Aminotab_productive = stats_IgG1_Aminotab[stats_IgG1_Aminotab[column_name] == "productive"]
            
            # Sort the filtered DataFrame by the selected column
            stats_IgG1_Aminotab_productive_sorted = stats_IgG1_Aminotab_productive.sort_values(by=column_name)
            
            
            # Create sets of IgE and IgG1 CDR3.IMGT for faster lookup
            stats_IgE_Aminotab_productive_CDR3IMGT = set(stats_IgE_Aminotab_productive['CDR3-IMGT'])

             # Check for the correct column name in stats_IgG1_Aminotab_productive
            if 'CDR3.IMGT' in stats_IgG1_Aminotab_productive.columns:
                column_cdr3 = 'CDR3.IMGT'
            elif 'CDR3-IMGT' in stats_IgG1_Aminotab_productive.columns:
                column_cdr3 = 'CDR3-IMGT'
            else:
                raise KeyError("Neither 'CDR3.IMGT' nor 'CDR3-IMGT' column is present in the DataFrame.")
            
            # Create a set of unique CDR3 sequences from the detected column
            stats_IgG1_Aminotab_productive_CDR3IMGT = set(stats_IgG1_Aminotab_productive[column_cdr3])

            
            # Create a copy of the data frame
            stats_tab_low_mutated_shared_IgE_VDJ_compare = stats_tab_low_mutated_shared_IgE_VDJ.copy()
            stats_tab_moderate_mutated_shared_IgE_VDJ_compare = stats_tab_moderate_mutated_shared_IgE_VDJ.copy()
            stats_tab_high_mutated_shared_IgE_VDJ_compare = stats_tab_high_mutated_shared_IgE_VDJ.copy()
            
            
            # Update the comparison IgE DataFrame with boolean masks
            stats_tab_low_mutated_shared_IgE_VDJ_compare['Match_in_low_mutated_IgG1_CDR3IMGT'] = stats_tab_low_mutated_shared_IgE_VDJ_compare['cdr3aa'].isin(stats_IgG1_Aminotab_productive_CDR3IMGT)
            stats_tab_moderate_mutated_shared_IgE_VDJ_compare['Match_in_moderate_mutated_IgG1_CDR3IMGT'] = stats_tab_moderate_mutated_shared_IgE_VDJ_compare['cdr3aa'].isin(stats_IgG1_Aminotab_productive_CDR3IMGT)
            stats_tab_high_mutated_shared_IgE_VDJ_compare['Match_in_high_mutated_IgG1_CDR3IMGT'] = stats_tab_high_mutated_shared_IgE_VDJ_compare['cdr3aa'].isin(stats_IgG1_Aminotab_productive_CDR3IMGT)
            
            
            # Ensure 'Match_in_low_moderate and high mutated_shared_IgG1_CDR3IMGT' is a string and remove leading/trailing spaces
            stats_tab_low_mutated_shared_IgE_VDJ_compare['Match_in_low_mutated_IgG1_CDR3IMGT'] = stats_tab_low_mutated_shared_IgE_VDJ_compare['Match_in_low_mutated_IgG1_CDR3IMGT'].astype(str).str.strip()
            stats_tab_moderate_mutated_shared_IgE_VDJ_compare['Match_in_moderate_mutated_IgG1_CDR3IMGT'] =stats_tab_moderate_mutated_shared_IgE_VDJ_compare['Match_in_moderate_mutated_IgG1_CDR3IMGT'].astype(str).str.strip()
            stats_tab_high_mutated_shared_IgE_VDJ_compare['Match_in_high_mutated_IgG1_CDR3IMGT'] = stats_tab_high_mutated_shared_IgE_VDJ_compare['Match_in_high_mutated_IgG1_CDR3IMGT'].astype(str).str.strip()
            
            
            # Now filter the divergent low_moderate and high mutated_shared_IgE clones
            divergent_tab_low_mutated_shared_IgE_VDJ = stats_tab_low_mutated_shared_IgE_VDJ_compare[(stats_tab_low_mutated_shared_IgE_VDJ_compare['Match_in_low_mutated_IgG1_CDR3IMGT'] == 'False')]
            divergent_tab_moderate_mutated_shared_IgE_VDJ = stats_tab_moderate_mutated_shared_IgE_VDJ_compare[(stats_tab_moderate_mutated_shared_IgE_VDJ_compare['Match_in_moderate_mutated_IgG1_CDR3IMGT'] == 'False')]
            divergent_tab_high_mutated_shared_IgE_VDJ = stats_tab_high_mutated_shared_IgE_VDJ_compare[(stats_tab_high_mutated_shared_IgE_VDJ_compare['Match_in_high_mutated_IgG1_CDR3IMGT'] == 'False')]
            
            
            # Now filter the non-divergent low_moderate and high mutated_shared_IgE clones
            non_divergent_tab_low_mutated_shared_IgE_VDJ = stats_tab_low_mutated_shared_IgE_VDJ_compare[(stats_tab_low_mutated_shared_IgE_VDJ_compare['Match_in_low_mutated_IgG1_CDR3IMGT'] == 'True')]
            non_divergent_tab_moderate_mutated_shared_IgE_VDJ = stats_tab_moderate_mutated_shared_IgE_VDJ_compare[(stats_tab_moderate_mutated_shared_IgE_VDJ_compare['Match_in_moderate_mutated_IgG1_CDR3IMGT'] == 'True')]
            non_divergent_tab_high_mutated_shared_IgE_VDJ = stats_tab_high_mutated_shared_IgE_VDJ_compare[(stats_tab_high_mutated_shared_IgE_VDJ_compare['Match_in_high_mutated_IgG1_CDR3IMGT'] == 'True')]
            
            
            #Filter the divergent low,moderate and high mutation with at least 2 CDR3 mutation
            # Create sets of highly, moderate and low mutated IgE seqids for faster lookup
            stats_IgE_low_mutated_seqid_CDR3_2 = set(divergent_tab_low_mutated_shared_IgE_VDJ['seqid'])
            stats_IgE_moderate_mutated_seqid_CDR3_2 = set(divergent_tab_moderate_mutated_shared_IgE_VDJ['seqid'])
            stats_IgE_high_mutated_seqid_CDR3_2 = set(divergent_tab_high_mutated_shared_IgE_VDJ['seqid'])
            
            
            # Create sets of highly, moderate and low mutated non divergent IgE seqids for faster lookup
            non_divergent_stats_IgE_low_mutated_seqid_CDR3_2 = set(non_divergent_tab_low_mutated_shared_IgE_VDJ['seqid'])
            non_divergent_stats_IgE_moderate_mutated_seqid_CDR3_2 = set(non_divergent_tab_moderate_mutated_shared_IgE_VDJ['seqid'])
            non_divergent_stats_IgE_high_mutated_seqid_CDR3_2 = set(non_divergent_tab_high_mutated_shared_IgE_VDJ['seqid'])
            
            
            # Create a copy of the data frame
            stats_IgE_mut_3_mut_level_compare = stats_IgE_mut_3.copy()
            
            
            # Update the comparison IgE DataFrame with boolean masks for the low,moderate and high divergent IgE subsets
            stats_IgE_mut_3_mut_level_compare['divergent_low_mutated_shared_IgE_seqid'] = stats_IgE_mut_3_mut_level_compare['Sequence ID'].isin(stats_IgE_low_mutated_seqid_CDR3_2)
            stats_IgE_mut_3_mut_level_compare['divergent_moderate_mutated_shared_IgE_seqid'] = stats_IgE_mut_3_mut_level_compare['Sequence ID'].isin(stats_IgE_moderate_mutated_seqid_CDR3_2)
            stats_IgE_mut_3_mut_level_compare['divergent_high_mutated_shared_IgE_seqid'] = stats_IgE_mut_3_mut_level_compare['Sequence ID'].isin(stats_IgE_high_mutated_seqid_CDR3_2)
            
            
            # Update the comparison IgE DataFrame with boolean masks for the low,moderate and high non-divergent IgE subsets
            stats_IgE_mut_3_mut_level_compare['non_divergent_low_mutated_shared_IgE_seqid'] = stats_IgE_mut_3_mut_level_compare['Sequence ID'].isin(non_divergent_stats_IgE_low_mutated_seqid_CDR3_2)
            stats_IgE_mut_3_mut_level_compare['non_divergent_moderate_mutated_shared_IgE_seqid'] = stats_IgE_mut_3_mut_level_compare['Sequence ID'].isin(non_divergent_stats_IgE_moderate_mutated_seqid_CDR3_2)
            stats_IgE_mut_3_mut_level_compare['non_divergent_high_mutated_shared_IgE_seqid'] = stats_IgE_mut_3_mut_level_compare['Sequence ID'].isin(non_divergent_stats_IgE_high_mutated_seqid_CDR3_2)
            
            
            # Ensure 'Match_in_low_moderate and high mutated_shared_divergent IgE_seqid' is a string and remove leading/trailing spaces
            stats_IgE_mut_3_mut_level_compare['divergent_low_mutated_shared_IgE_seqid'] = stats_IgE_mut_3_mut_level_compare['divergent_low_mutated_shared_IgE_seqid'].astype(str).str.strip()
            stats_IgE_mut_3_mut_level_compare['divergent_moderate_mutated_shared_IgE_seqid'] =stats_IgE_mut_3_mut_level_compare['divergent_moderate_mutated_shared_IgE_seqid'].astype(str).str.strip()
            stats_IgE_mut_3_mut_level_compare['divergent_high_mutated_shared_IgE_seqid'] = stats_IgE_mut_3_mut_level_compare['divergent_high_mutated_shared_IgE_seqid'].astype(str).str.strip()
            
            
            # Ensure 'Match_in_low_moderate and high mutated_non divergent shared_IgE_seqid' is a string and remove leading/trailing spaces
            stats_IgE_mut_3_mut_level_compare['non_divergent_low_mutated_shared_IgE_seqid'] = stats_IgE_mut_3_mut_level_compare['non_divergent_low_mutated_shared_IgE_seqid'].astype(str).str.strip()
            stats_IgE_mut_3_mut_level_compare['non_divergent_moderate_mutated_shared_IgE_seqid'] =stats_IgE_mut_3_mut_level_compare['non_divergent_moderate_mutated_shared_IgE_seqid'].astype(str).str.strip()
            stats_IgE_mut_3_mut_level_compare['non_divergent_high_mutated_shared_IgE_seqid'] = stats_IgE_mut_3_mut_level_compare['non_divergent_high_mutated_shared_IgE_seqid'].astype(str).str.strip()
            
            
            # Now filter the divergent low_moderate and high mutated_shared_IgE mutation table with  > 2 mutations
            divergent_tab_low_mutated_shared_IgE_VDJ_CDR3_2 = stats_IgE_mut_3_mut_level_compare[(stats_IgE_mut_3_mut_level_compare['divergent_low_mutated_shared_IgE_seqid'] == 'True') & (stats_IgE_mut_3_mut_level_compare['CDR3-IMGT Nb of nonsilent mutations-ex'] >= 2)]
            divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3_2 = stats_IgE_mut_3_mut_level_compare[(stats_IgE_mut_3_mut_level_compare['divergent_moderate_mutated_shared_IgE_seqid'] == 'True') & (stats_IgE_mut_3_mut_level_compare['CDR3-IMGT Nb of nonsilent mutations-ex'] >= 2)]
            divergent_tab_high_mutated_shared_IgE_VDJ_CDR3_2 = stats_IgE_mut_3_mut_level_compare[(stats_IgE_mut_3_mut_level_compare['divergent_high_mutated_shared_IgE_seqid'] == 'True')& (stats_IgE_mut_3_mut_level_compare['CDR3-IMGT Nb of nonsilent mutations-ex'] >= 2)]
            
            # Now filter the non-divergent low_moderate and high mutated_shared_IgE mutation table with  > 2 mutations
            non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3_2 = stats_IgE_mut_3_mut_level_compare[(stats_IgE_mut_3_mut_level_compare['non_divergent_low_mutated_shared_IgE_seqid'] == 'True') & (stats_IgE_mut_3_mut_level_compare['CDR3-IMGT Nb of nonsilent mutations-ex'] >= 2)]
            non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3_2 = stats_IgE_mut_3_mut_level_compare[(stats_IgE_mut_3_mut_level_compare['non_divergent_moderate_mutated_shared_IgE_seqid'] == 'True') & (stats_IgE_mut_3_mut_level_compare['CDR3-IMGT Nb of nonsilent mutations-ex'] >= 2)]
            non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3_2 = stats_IgE_mut_3_mut_level_compare[(stats_IgE_mut_3_mut_level_compare['non_divergent_high_mutated_shared_IgE_seqid'] == 'True')& (stats_IgE_mut_3_mut_level_compare['CDR3-IMGT Nb of nonsilent mutations-ex'] >= 2)]
            
            
            # Create sets of IgE CDR3.2.seqid for faster lookup for the divergent IgEs
            divergent_tab_low_mutated_shared_IgE_VDJ_CDR3_2_seqid = set(divergent_tab_low_mutated_shared_IgE_VDJ_CDR3_2['Sequence ID'])
            divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3_2_seqid = set(divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3_2['Sequence ID'])
            divergent_tab_high_mutated_shared_IgE_VDJ_CDR3_2_seqid = set(divergent_tab_high_mutated_shared_IgE_VDJ_CDR3_2['Sequence ID'])
            
            # Create sets of IgE CDR3.2.seqid for faster lookup for the non-divergent IgEs
            non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3_2_seqid = set(non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3_2['Sequence ID'])
            non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3_2_seqid = set(non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3_2['Sequence ID'])
            non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3_2_seqid = set(non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3_2['Sequence ID'])
            
            # Create a copy of the data frame for the divergent IgEs
            divergent_tab_low_mutated_shared_IgE_VDJ_compare = divergent_tab_low_mutated_shared_IgE_VDJ.copy()
            divergent_tab_moderate_mutated_shared_IgE_VDJ_compare = divergent_tab_moderate_mutated_shared_IgE_VDJ.copy()
            divergent_tab_high_mutated_shared_IgE_VDJ_compare = divergent_tab_high_mutated_shared_IgE_VDJ.copy()
            
            # Create a copy of the data frame for the non-divergent IgEs
            non_divergent_tab_low_mutated_shared_IgE_VDJ_compare = non_divergent_tab_low_mutated_shared_IgE_VDJ.copy()
            non_divergent_tab_moderate_mutated_shared_IgE_VDJ_compare = non_divergent_tab_moderate_mutated_shared_IgE_VDJ.copy()
            non_divergent_tab_high_mutated_shared_IgE_VDJ_compare = non_divergent_tab_high_mutated_shared_IgE_VDJ.copy()
            
            # Update the comparison IgE DataFrame with boolean masks for the divergent IgEs
            divergent_tab_low_mutated_shared_IgE_VDJ_compare['Match_in_low_mutated_shared_IgE_VDJ_CDR3_2_seqid'] = divergent_tab_low_mutated_shared_IgE_VDJ_compare['seqid'].isin(divergent_tab_low_mutated_shared_IgE_VDJ_CDR3_2_seqid)
            divergent_tab_moderate_mutated_shared_IgE_VDJ_compare['Match_in_moderate_mutated_shared_IgE_VDJ_CDR3_2_seqid'] = divergent_tab_moderate_mutated_shared_IgE_VDJ_compare['seqid'].isin(divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3_2_seqid)
            divergent_tab_high_mutated_shared_IgE_VDJ_compare['Match_in_high_mutated_shared_IgE_VDJ_CDR3_2_seqid'] = divergent_tab_high_mutated_shared_IgE_VDJ_compare['seqid'].isin(divergent_tab_high_mutated_shared_IgE_VDJ_CDR3_2_seqid)
            
            # Update the comparison IgE DataFrame with boolean masks for the non-divergent IgEs
            non_divergent_tab_low_mutated_shared_IgE_VDJ_compare['Match_in_low_mutated_shared_IgE_VDJ_CDR3_2_seqid'] = non_divergent_tab_low_mutated_shared_IgE_VDJ_compare['seqid'].isin(non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3_2_seqid)
            non_divergent_tab_moderate_mutated_shared_IgE_VDJ_compare['Match_in_moderate_mutated_shared_IgE_VDJ_CDR3_2_seqid'] = non_divergent_tab_moderate_mutated_shared_IgE_VDJ_compare['seqid'].isin(non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3_2_seqid)
            non_divergent_tab_high_mutated_shared_IgE_VDJ_compare['Match_in_high_mutated_shared_IgE_VDJ_CDR3_2_seqid'] = non_divergent_tab_high_mutated_shared_IgE_VDJ_compare['seqid'].isin(non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3_2_seqid)
            
            # Ensure 'Match_in_the divergent low_moderate and high mutated_shared_IgE_seqid' is a string and remove leading/trailing spaces
            divergent_tab_low_mutated_shared_IgE_VDJ_compare['Match_in_low_mutated_shared_IgE_VDJ_CDR3_2_seqid'] = divergent_tab_low_mutated_shared_IgE_VDJ_compare['Match_in_low_mutated_shared_IgE_VDJ_CDR3_2_seqid'].astype(str).str.strip()
            divergent_tab_moderate_mutated_shared_IgE_VDJ_compare['Match_in_moderate_mutated_shared_IgE_VDJ_CDR3_2_seqid'] = divergent_tab_moderate_mutated_shared_IgE_VDJ_compare['Match_in_moderate_mutated_shared_IgE_VDJ_CDR3_2_seqid'].astype(str).str.strip()
            divergent_tab_high_mutated_shared_IgE_VDJ_compare['Match_in_high_mutated_shared_IgE_VDJ_CDR3_2_seqid'] = divergent_tab_high_mutated_shared_IgE_VDJ_compare['Match_in_high_mutated_shared_IgE_VDJ_CDR3_2_seqid'].astype(str).str.strip()
            
            # Ensure 'Match_in_the non-divergent low_moderate and high mutated_shared_IgE_seqid' is a string and remove leading/trailing spaces
            non_divergent_tab_low_mutated_shared_IgE_VDJ_compare['Match_in_low_mutated_shared_IgE_VDJ_CDR3_2_seqid'] = non_divergent_tab_low_mutated_shared_IgE_VDJ_compare['Match_in_low_mutated_shared_IgE_VDJ_CDR3_2_seqid'].astype(str).str.strip()
            non_divergent_tab_moderate_mutated_shared_IgE_VDJ_compare['Match_in_moderate_mutated_shared_IgE_VDJ_CDR3_2_seqid'] = non_divergent_tab_moderate_mutated_shared_IgE_VDJ_compare['Match_in_moderate_mutated_shared_IgE_VDJ_CDR3_2_seqid'].astype(str).str.strip()
            non_divergent_tab_high_mutated_shared_IgE_VDJ_compare['Match_in_high_mutated_shared_IgE_VDJ_CDR3_2_seqid'] = non_divergent_tab_high_mutated_shared_IgE_VDJ_compare['Match_in_high_mutated_shared_IgE_VDJ_CDR3_2_seqid'].astype(str).str.strip()
            
            # Now filter the divergent low_moderate and high mutated_shared_IgE clones
            divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone = divergent_tab_low_mutated_shared_IgE_VDJ_compare[(divergent_tab_low_mutated_shared_IgE_VDJ_compare['Match_in_low_mutated_shared_IgE_VDJ_CDR3_2_seqid'] == 'True')]
            divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone = divergent_tab_moderate_mutated_shared_IgE_VDJ_compare[(divergent_tab_moderate_mutated_shared_IgE_VDJ_compare['Match_in_moderate_mutated_shared_IgE_VDJ_CDR3_2_seqid'] == 'True')]
            divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone = divergent_tab_high_mutated_shared_IgE_VDJ_compare[(divergent_tab_high_mutated_shared_IgE_VDJ_compare['Match_in_high_mutated_shared_IgE_VDJ_CDR3_2_seqid'] == 'True')]
            
            
            # Now filter the non-divergent low_moderate and high mutated_shared_IgE clones
            non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone = non_divergent_tab_low_mutated_shared_IgE_VDJ_compare[(non_divergent_tab_low_mutated_shared_IgE_VDJ_compare['Match_in_low_mutated_shared_IgE_VDJ_CDR3_2_seqid'] == 'True')]
            non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone = non_divergent_tab_moderate_mutated_shared_IgE_VDJ_compare[(non_divergent_tab_moderate_mutated_shared_IgE_VDJ_compare['Match_in_moderate_mutated_shared_IgE_VDJ_CDR3_2_seqid'] == 'True')]
            non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone = non_divergent_tab_high_mutated_shared_IgE_VDJ_compare[(non_divergent_tab_high_mutated_shared_IgE_VDJ_compare['Match_in_high_mutated_shared_IgE_VDJ_CDR3_2_seqid'] == 'True')]
            
            
            # Determine the number of divergent low, moderate, and high mutated IgE clones
            unique_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone = divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone['IgE_VDJ'].nunique(dropna=True)
            print(f'Number of divergent CDR3 low mutated IgE clones: {unique_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone}')
            
            unique_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone = divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone['IgE_VDJ'].nunique(dropna=True)
            print(f'Number of divergent CDR3 moderate mutated IgE clones: {unique_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone}')
            
            unique_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone = divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone['IgE_VDJ'].nunique(dropna=True)
            print(f'Number of divergent CDR3 high mutated IgE clones: {unique_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone}')
            
            
            # Determine the number of non-divergent low, moderate, and high mutated IgE clones
            unique_non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone = non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone['IgE_VDJ'].nunique(dropna=True)
            print(f'Number of non_divergent CDR3 low mutated IgE clones: {unique_non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone}')
            
            unique_non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone = non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone['IgE_VDJ'].nunique(dropna=True)
            print(f'Number of non_divergent CDR3 moderate mutated IgE clones: {unique_non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone}')
            
            unique_non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone = non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone['IgE_VDJ'].nunique(dropna=True)
            print(f'Number of non_divergent CDR3 high mutated IgE clones: {unique_non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone}')
            
            
            
            # Determine the average copy number of divergent and non-divergent low, moderate and high mutated IgE clones
            
            # Calculate the mean copynumber for the divergent low mutated IgE clones, ensuring NaN results are handled
            mean_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone = divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone['onecopy'].mean()
            mean_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone = 0 if pd.isna(mean_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone) else mean_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone
            print(f'mean copynumber of divergent CDR3 low mutated IgE: {mean_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone}')
            
            # Calculate the mean copynumber for the divergent moderate mutated IgE clones, ensuring NaN results are handled
            mean_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone = divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone['onecopy'].mean()
            mean_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone = 0 if pd.isna(mean_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone) else mean_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone
            print(f'mean copynumber of divergent CDR3 moderate mutated IgE: {mean_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone}')
            
            # Calculate the mean copynumber for the divergent high mutated IgE clones, ensuring NaN results are handled
            mean_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone = divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone['onecopy'].mean()
            mean_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone = 0 if pd.isna(mean_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone) else mean_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone
            print(f'mean copynumber of divergent CDR3 high mutated IgE: {mean_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone}')
            
            
            # Calculate the mean copynumber for non_divergent low mutated IgE clones, ensuring NaN results are handled
            mean_non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone = non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone['onecopy'].mean()
            mean_non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone = 0 if pd.isna(mean_non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone) else mean_non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone
            print(f'mean copynumber of non_divergent CDR3 low mutated IgE: {mean_non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone}')
            
            # Calculate the mean copynumber for non_divergent moderate mutated IgE clones, ensuring NaN results are handled
            mean_non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone = non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone['onecopy'].mean()
            mean_non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone = 0 if pd.isna(mean_non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone) else mean_non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone
            print(f'mean copynumber of non_divergent CDR3 moderate mutated IgE: {mean_non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone}')
            
            # Calculate the mean copynumber for non_divergent high mutated IgE clones, ensuring NaN results are handled
            mean_non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone = non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone['onecopy'].mean()
            mean_non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone = 0 if pd.isna(mean_non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone) else mean_non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone
            print(f'mean copynumber of non_divergent CDR3 high mutated IgE: {mean_non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone}')
            
            
            # Calculate the sum copynumber for low mutated IgE clones, ensuring NaN results are handled
            sum_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone = divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone['onecopy'].sum()
            sum_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone = 0 if pd.isna(sum_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone) else sum_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone
            print(f'sum copynumber of divergent CDR3 low mutated IgE: {sum_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone}')
            
            # Calculate the mean copynumber for moderate mutated IgE clones, ensuring NaN results are handled
            sum_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone = divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone['onecopy'].sum()
            sum_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone = 0 if pd.isna(sum_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone) else sum_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone
            print(f'sum copynumber of divergent CDR3 moderate mutated IgE: {sum_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone}')
            
            # Calculate the mean copynumber for high mutated IgE clones, ensuring NaN results are handled
            sum_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone = divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone['onecopy'].sum()
            sum_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone = 0 if pd.isna(sum_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone) else sum_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone
            print(f'sum copynumber of divergent CDR3 high mutated IgE: {sum_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone}')
            
            
            
            # Calculate the sum copynumber for low mutated IgE clones, ensuring NaN results are handled
            sum_non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone = non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone['onecopy'].sum()
            sum_non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone = 0 if pd.isna(sum_non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone) else sum_non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone
            print(f'sum copynumber of non_divergent CDR3 low mutated IgE: {sum_non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone}')
            
            # Calculate the mean copynumber for moderate mutated IgE clones, ensuring NaN results are handled
            sum_non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone = non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone['onecopy'].sum()
            sum_non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone = 0 if pd.isna(sum_non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone) else sum_non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone
            print(f'sum copynumber of non_divergent CDR3 moderate mutated IgE: {sum_non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone}')
            
            # Calculate the mean copynumber for high mutated IgE clones, ensuring NaN results are handled
            sum_non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone = non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone['onecopy'].sum()
            sum_non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone = 0 if pd.isna(sum_non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone) else sum_non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone
            print(f'sum copynumber of non_divergent CDR3 high mutated IgE: {sum_non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone}')
            
            
            
            
            # Filter out the clonally related divergent low, moderate and high mutated IgG1 counterpart
            # Create sets of low, moderate and high IgE CDR3.2.seqid for faster lookup
            divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone_seqid = set(divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone['IgE_VDJ'])
            divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone_seqid = set(divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone['IgE_VDJ'])
            divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone_seqid = set(divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone['IgE_VDJ'])
            
            
            # Create a copy of the data frame
            clonally_related_low_mod_high_shared_IgG1 = shared_IgG1.copy()
            
            
            # Update the comparison IgG1 DataFrame with boolean masks
            clonally_related_low_mod_high_shared_IgG1['divergent_low_mutated_shared_IgE_VDJ_in_IgG1_seqid'] = clonally_related_low_mod_high_shared_IgG1['IgG1_VDJ'].isin(divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone_seqid)
            clonally_related_low_mod_high_shared_IgG1['divergent_moderate_mutated_shared_IgE_VDJ_in_IgG1_seqid'] = clonally_related_low_mod_high_shared_IgG1['IgG1_VDJ'].isin(divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone_seqid)
            clonally_related_low_mod_high_shared_IgG1['divergent_high_mutated_shared_IgE_VDJ_in_IgG1_seqid'] = clonally_related_low_mod_high_shared_IgG1['IgG1_VDJ'].isin(divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone_seqid)
            
            # Ensure 'Match_in_low_moderate and high mutated_shared_IgE_seqid' is a string and remove leading/trailing spaces
            clonally_related_low_mod_high_shared_IgG1['divergent_low_mutated_shared_IgE_VDJ_in_IgG1_seqid'] = clonally_related_low_mod_high_shared_IgG1['divergent_low_mutated_shared_IgE_VDJ_in_IgG1_seqid'].astype(str).str.strip()
            clonally_related_low_mod_high_shared_IgG1['divergent_moderate_mutated_shared_IgE_VDJ_in_IgG1_seqid'] =clonally_related_low_mod_high_shared_IgG1['divergent_moderate_mutated_shared_IgE_VDJ_in_IgG1_seqid'].astype(str).str.strip()
            clonally_related_low_mod_high_shared_IgG1['divergent_high_mutated_shared_IgE_VDJ_in_IgG1_seqid'] = clonally_related_low_mod_high_shared_IgG1['divergent_high_mutated_shared_IgE_VDJ_in_IgG1_seqid'].astype(str).str.strip()
            
            
            # Now filter the low_moderate and high mutated_shared_clonally related IgG1 clones
            clonally_related_low_mutated_shared_IgG1 = clonally_related_low_mod_high_shared_IgG1[(clonally_related_low_mod_high_shared_IgG1['divergent_low_mutated_shared_IgE_VDJ_in_IgG1_seqid'] == 'True')]
            clonally_related_moderate_mutated_shared_IgG1 = clonally_related_low_mod_high_shared_IgG1[(clonally_related_low_mod_high_shared_IgG1['divergent_moderate_mutated_shared_IgE_VDJ_in_IgG1_seqid'] == 'True')]
            clonally_related_high_mutated_shared_IgG1 = clonally_related_low_mod_high_shared_IgG1[(clonally_related_low_mod_high_shared_IgG1['divergent_high_mutated_shared_IgE_VDJ_in_IgG1_seqid'] == 'True')]
            
            # Determine the average copy number of divergent low, moderate and high mutated clonally related IgG1 counterpart
            
            # Calculate the mean copynumber for low mutated IgE clones IgG1 counterpart, ensuring NaN results are handled
            mean_clonally_related_low_mutated_shared_IgG1 = clonally_related_low_mutated_shared_IgG1['onecopy'].mean()
            mean_clonally_related_low_mutated_shared_IgG1 = 0 if pd.isna(mean_clonally_related_low_mutated_shared_IgG1) else mean_clonally_related_low_mutated_shared_IgG1
            print(f'mean copynumber of IgG1 clonally related to low mutated IgE: {mean_clonally_related_low_mutated_shared_IgG1}')
            
            # Calculate the mean copynumber for moderate mutated IgE clones IgG1 counterpart, ensuring NaN results are handled
            mean_clonally_related_moderate_mutated_shared_IgG1 = clonally_related_moderate_mutated_shared_IgG1['onecopy'].mean()
            mean_clonally_related_moderate_mutated_shared_IgG1 = 0 if pd.isna(mean_clonally_related_moderate_mutated_shared_IgG1) else mean_clonally_related_moderate_mutated_shared_IgG1
            print(f'mean copynumber of IgG1 clonally related to moderate mutated IgE: {mean_clonally_related_moderate_mutated_shared_IgG1}')
            
            # Calculate the mean copynumber for moderate mutated IgE clones IgG1 counterpart, ensuring NaN results are handled
            mean_clonally_related_high_mutated_shared_IgG1 = clonally_related_high_mutated_shared_IgG1['onecopy'].mean()
            mean_clonally_related_high_mutated_shared_IgG1 = 0 if pd.isna(mean_clonally_related_high_mutated_shared_IgG1) else mean_clonally_related_high_mutated_shared_IgG1
            print(f'mean copynumber of IgG1 clonally related to high mutated IgE: {mean_clonally_related_high_mutated_shared_IgG1}')
            
            
            # Calculate the sum copynumber for low mutated IgE clones IgG1 counterpart, ensuring NaN results are handled
            sum_clonally_related_low_mutated_shared_IgG1 = clonally_related_low_mutated_shared_IgG1['onecopy'].sum()
            sum_clonally_related_low_mutated_shared_IgG1 = 0 if pd.isna(sum_clonally_related_low_mutated_shared_IgG1) else sum_clonally_related_low_mutated_shared_IgG1
            print(f'sum copynumber of IgG1 clonally related to low mutated IgE: {sum_clonally_related_low_mutated_shared_IgG1}')
            
            # Calculate the mean copynumber for moderate mutated IgE clones IgG1 counterpart, ensuring NaN results are handled
            sum_clonally_related_moderate_mutated_shared_IgG1 = clonally_related_moderate_mutated_shared_IgG1['onecopy'].sum()
            sum_clonally_related_moderate_mutated_shared_IgG1 = 0 if pd.isna(sum_clonally_related_moderate_mutated_shared_IgG1) else sum_clonally_related_moderate_mutated_shared_IgG1
            print(f'sum copynumber of IgG1 clonally related to moderate mutated IgE: {sum_clonally_related_moderate_mutated_shared_IgG1}')
            
            # Calculate the mean copynumber for moderate mutated IgE clones IgG1 counterpart, ensuring NaN results are handled
            sum_clonally_related_high_mutated_shared_IgG1 = clonally_related_high_mutated_shared_IgG1['onecopy'].sum()
            sum_clonally_related_high_mutated_shared_IgG1 = 0 if pd.isna(sum_clonally_related_high_mutated_shared_IgG1) else sum_clonally_related_high_mutated_shared_IgG1
            print(f'sum copynumber of IgG1 clonally related to high mutated IgE: {sum_clonally_related_high_mutated_shared_IgG1}')
            
            
            
            # Analyse the low, moderate and highly divergent based on the CDR3-CDR2 aminoacid sequence
            # Merge the CDR3 and CDR2 aminoacid sequence of Total IgG1
            stats_IgE_Aminotab_productive_CDR3CDR2 = stats_IgE_Aminotab_productive.copy()
            stats_IgE_Aminotab_productive_CDR3CDR2 ['CDR3CDR2-IMGT'] = stats_IgE_Aminotab_productive_CDR3CDR2['CDR2-IMGT'].astype(str) + stats_IgE_Aminotab_productive_CDR3CDR2['CDR3-IMGT'].astype(str)
            
            # Remove string spaces in each row
            stats_IgE_Aminotab_productive_CDR3CDR2 ['CDR3CDR2-IMGT'] = stats_IgE_Aminotab_productive_CDR3CDR2 ['CDR3CDR2-IMGT'].str.replace(' ', '', regex=False)
            
            # Merge the CDR3 and CDR2 aminoacid sequence of Total IgG1
            stats_IgG1_Aminotab_productive_CDR3CDR2 = stats_IgG1_Aminotab_productive.copy()
            
                        # Check for the correct column names in stats_IgG1_Aminotab_productive_CDR3CDR2
            if 'CDR2.IMGT' in stats_IgG1_Aminotab_productive_CDR3CDR2.columns:
                column_cdr2 = 'CDR2.IMGT'
            elif 'CDR2-IMGT' in stats_IgG1_Aminotab_productive_CDR3CDR2.columns:
                column_cdr2 = 'CDR2-IMGT'
            else:
                raise KeyError("Neither 'CDR2.IMGT' nor 'CDR2-IMGT' column is present in the DataFrame.")
            
            if 'CDR3.IMGT' in stats_IgG1_Aminotab_productive_CDR3CDR2.columns:
                column_cdr3 = 'CDR3.IMGT'
            elif 'CDR3-IMGT' in stats_IgG1_Aminotab_productive_CDR3CDR2.columns:
                column_cdr3 = 'CDR3-IMGT'
            else:
                raise KeyError("Neither 'CDR3.IMGT' nor 'CDR3-IMGT' column is present in the DataFrame.")
            
            # Concatenate the CDR2 and CDR3 columns as strings
            stats_IgG1_Aminotab_productive_CDR3CDR2['CDR3CDR2-IMGT'] = \
                stats_IgG1_Aminotab_productive_CDR3CDR2[column_cdr2].astype(str) + stats_IgG1_Aminotab_productive_CDR3CDR2[column_cdr3].astype(str)

            
            # Remove string spaces in each row
            stats_IgG1_Aminotab_productive_CDR3CDR2 ['CDR3CDR2-IMGT'] = stats_IgG1_Aminotab_productive_CDR3CDR2 ['CDR3CDR2-IMGT'].str.replace(' ', '', regex=False)
            
            
            #Extract the high, moderate and low mutated CDR3CDR2
            
            # Update the comparison IgE DataFrame with boolean masks
            stats_IgE_Aminotab_productive_CDR3CDR2['low_mutated_CDR3CDR2-IMGT'] = stats_IgE_Aminotab_productive_CDR3CDR2['Sequence ID'].isin(low_mutated_seqids_Shared_IgE)
            stats_IgE_Aminotab_productive_CDR3CDR2['moderate_mutated_CDR3CDR2-IMGT'] = stats_IgE_Aminotab_productive_CDR3CDR2['Sequence ID'].isin(moderate_mutated_seqids_Shared_IgE)
            stats_IgE_Aminotab_productive_CDR3CDR2['high_mutated_CDR3CDR2-IMGT'] = stats_IgE_Aminotab_productive_CDR3CDR2['Sequence ID'].isin(High_mutated_seqids_Shared_IgE)
            
            # Ensure 'Match_in_low_moderate and high mutated_shared_IgE_seqid' is a string and remove leading/trailing spaces
            stats_IgE_Aminotab_productive_CDR3CDR2 ['low_mutated_CDR3CDR2-IMGT'] = stats_IgE_Aminotab_productive_CDR3CDR2 ['low_mutated_CDR3CDR2-IMGT'].astype(str).str.strip()
            stats_IgE_Aminotab_productive_CDR3CDR2['moderate_mutated_CDR3CDR2-IMGT'] =stats_IgE_Aminotab_productive_CDR3CDR2 ['moderate_mutated_CDR3CDR2-IMGT'].astype(str).str.strip()
            stats_IgE_Aminotab_productive_CDR3CDR2 ['high_mutated_CDR3CDR2-IMGT'] = stats_IgE_Aminotab_productive_CDR3CDR2 ['high_mutated_CDR3CDR2-IMGT'].astype(str).str.strip()
            
            # Now filter the low_moderate and high mutated_shared_IgE CDR3CDR2
            stats_low_mutated_IgE_Aminotab_productive_CDR3CDR2 = stats_IgE_Aminotab_productive_CDR3CDR2[(stats_IgE_Aminotab_productive_CDR3CDR2['low_mutated_CDR3CDR2-IMGT'] == 'True')]
            stats_moderate_mutated_IgE_Aminotab_productive_CDR3CDR2 = stats_IgE_Aminotab_productive_CDR3CDR2[(stats_IgE_Aminotab_productive_CDR3CDR2['moderate_mutated_CDR3CDR2-IMGT'] == 'True')]
            stats_high_mutated_IgE_Aminotab_productive_CDR3CDR2 = stats_IgE_Aminotab_productive_CDR3CDR2[(stats_IgE_Aminotab_productive_CDR3CDR2['high_mutated_CDR3CDR2-IMGT'] == 'True')]
            
            #Compare the IgG1 CDR3CDR2 aminoacid table to the high, moderate and high mutated CDR3CDR2
            # Create sets of IgG1 CDR3CDR2 for faster lookup
            stats_IgG1_Aminotab_productive_CDR3CDR2_IMGTAA = set(stats_IgG1_Aminotab_productive_CDR3CDR2['CDR3CDR2-IMGT'])
            
            # Ensure explicit copies of sliced DataFrames
            stats_low_mutated_IgE_Aminotab_productive_CDR3CDR2 = stats_low_mutated_IgE_Aminotab_productive_CDR3CDR2.copy()
            stats_moderate_mutated_IgE_Aminotab_productive_CDR3CDR2 = stats_moderate_mutated_IgE_Aminotab_productive_CDR3CDR2.copy()
            stats_high_mutated_IgE_Aminotab_productive_CDR3CDR2 = stats_high_mutated_IgE_Aminotab_productive_CDR3CDR2.copy()
            
            # Modify the columns in the copied DataFrames
            stats_low_mutated_IgE_Aminotab_productive_CDR3CDR2['low_mutated_CDR3CDR2-in_IgG1_CDR3CDR2'] = stats_low_mutated_IgE_Aminotab_productive_CDR3CDR2['CDR3CDR2-IMGT'].isin(stats_IgG1_Aminotab_productive_CDR3CDR2_IMGTAA)
            stats_moderate_mutated_IgE_Aminotab_productive_CDR3CDR2['moderate_mutated_CDR3CDR2-in_IgG1_CDR3CDR2'] = stats_moderate_mutated_IgE_Aminotab_productive_CDR3CDR2['CDR3CDR2-IMGT'].isin(stats_IgG1_Aminotab_productive_CDR3CDR2_IMGTAA)
            stats_high_mutated_IgE_Aminotab_productive_CDR3CDR2['high_mutated_CDR3CDR2-in_IgG1_CDR3CDR2'] = stats_high_mutated_IgE_Aminotab_productive_CDR3CDR2['CDR3CDR2-IMGT'].isin(stats_IgG1_Aminotab_productive_CDR3CDR2_IMGTAA)
            
            # Filter the divergent low, moderate, and high mutated shared IgE CDR3CDR2
            divergent_stats_low_mutated_IgE_Aminotab_productive_CDR3CDR2 = stats_low_mutated_IgE_Aminotab_productive_CDR3CDR2[~stats_low_mutated_IgE_Aminotab_productive_CDR3CDR2['low_mutated_CDR3CDR2-in_IgG1_CDR3CDR2']]
            divergent_stats_moderate_mutated_IgE_Aminotab_productive_CDR3CDR2 = stats_moderate_mutated_IgE_Aminotab_productive_CDR3CDR2[~stats_moderate_mutated_IgE_Aminotab_productive_CDR3CDR2['moderate_mutated_CDR3CDR2-in_IgG1_CDR3CDR2']]
            divergent_stats_high_mutated_IgE_Aminotab_productive_CDR3CDR2 = stats_high_mutated_IgE_Aminotab_productive_CDR3CDR2[~stats_high_mutated_IgE_Aminotab_productive_CDR3CDR2['high_mutated_CDR3CDR2-in_IgG1_CDR3CDR2']]
            
            
            # Extract the divergent low_moderate and high mutated_shared_IgE CDR3CDR2 mutation clones with >= 2 CDR3 non silent mutation
            # Create sets of IgE and IgG1 CDR3.IMGT for faster lookup
            divergent_stats_low_mutated_IgE_Aminotab_productive_CDR3CDR2_seqid = set(divergent_stats_low_mutated_IgE_Aminotab_productive_CDR3CDR2['Sequence ID'])
            divergent_stats_moderate_mutated_IgE_Aminotab_productive_CDR3CDR2_seqid = set(divergent_stats_moderate_mutated_IgE_Aminotab_productive_CDR3CDR2['Sequence ID'])
            divergent_stats_high_mutated_IgE_Aminotab_productive_CDR3CDR2_seqid = set(divergent_stats_high_mutated_IgE_Aminotab_productive_CDR3CDR2['Sequence ID'])
            
            # Create a copy of the data frame
            stats_IgE_mut_3_CDR3CDR2_CDR3_2 = stats_IgE_mut_3.copy()
            
            # Update the comparison IgE DataFrame with boolean masks
            stats_IgE_mut_3_CDR3CDR2_CDR3_2['Match_in_divergent_low_mutated_IgE_seqid'] = stats_IgE_mut_3_CDR3CDR2_CDR3_2['Sequence ID'].isin(divergent_stats_low_mutated_IgE_Aminotab_productive_CDR3CDR2_seqid)
            stats_IgE_mut_3_CDR3CDR2_CDR3_2['Match_in_divergent_moderate_mutated_IgE_seqid'] = stats_IgE_mut_3_CDR3CDR2_CDR3_2['Sequence ID'].isin(divergent_stats_moderate_mutated_IgE_Aminotab_productive_CDR3CDR2_seqid)
            stats_IgE_mut_3_CDR3CDR2_CDR3_2['Match_in_divergent_high_mutated_IgE_seqid'] = stats_IgE_mut_3_CDR3CDR2_CDR3_2['Sequence ID'].isin(divergent_stats_high_mutated_IgE_Aminotab_productive_CDR3CDR2_seqid)
            
            # Ensure 'Match_in_low_moderate and high mutated_shared_IgG1_CDR3IMGT' is a string and remove leading/trailing spaces
            stats_IgE_mut_3_CDR3CDR2_CDR3_2['Match_in_divergent_low_mutated_IgE_seqid'] = stats_IgE_mut_3_CDR3CDR2_CDR3_2['Match_in_divergent_low_mutated_IgE_seqid'].astype(str).str.strip()
            stats_IgE_mut_3_CDR3CDR2_CDR3_2['Match_in_divergent_moderate_mutated_IgE_seqid'] =stats_IgE_mut_3_CDR3CDR2_CDR3_2['Match_in_divergent_moderate_mutated_IgE_seqid'].astype(str).str.strip()
            stats_IgE_mut_3_CDR3CDR2_CDR3_2['Match_in_divergent_high_mutated_IgE_seqid'] = stats_IgE_mut_3_CDR3CDR2_CDR3_2['Match_in_divergent_high_mutated_IgE_seqid'].astype(str).str.strip()
            
            # Now filter the low_moderate and high mutated_shared_IgE mutation table
            divergent_tab_low_mutated_shared_IgE_CDR3CDR2_mutab = stats_IgE_mut_3_CDR3CDR2_CDR3_2[(stats_IgE_mut_3_CDR3CDR2_CDR3_2['Match_in_divergent_low_mutated_IgE_seqid'] == 'True') & (stats_IgE_mut_3_CDR3CDR2_CDR3_2['CDR3-IMGT Nb of nonsilent mutations-ex'] >= 2)]
            divergent_tab_moderate_mutated_shared_IgE_CDR3CDR2_mutab = stats_IgE_mut_3_CDR3CDR2_CDR3_2[(stats_IgE_mut_3_CDR3CDR2_CDR3_2['Match_in_divergent_moderate_mutated_IgE_seqid'] == 'True')& (stats_IgE_mut_3_CDR3CDR2_CDR3_2['CDR3-IMGT Nb of nonsilent mutations-ex'] >= 2)]
            divergent_tab_high_mutated_shared_IgE_CDR3CDR2_mutab = stats_IgE_mut_3_CDR3CDR2_CDR3_2[(stats_IgE_mut_3_CDR3CDR2_CDR3_2['Match_in_divergent_high_mutated_IgE_seqid'] == 'True')& (stats_IgE_mut_3_CDR3CDR2_CDR3_2['CDR3-IMGT Nb of nonsilent mutations-ex'] >= 2)]
            
            # Create sets of IgE CDR3.2.seqid for faster lookup
            divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2_2_seqid = set(divergent_tab_low_mutated_shared_IgE_CDR3CDR2_mutab['Sequence ID'])
            divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2_2_seqid = set(divergent_tab_moderate_mutated_shared_IgE_CDR3CDR2_mutab['Sequence ID'])
            divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2_2_seqid = set(divergent_tab_high_mutated_shared_IgE_CDR3CDR2_mutab['Sequence ID'])
            
            # Create a copy of the data frame
            divergent_tab_low_moderate_high_mutated_shared_IgE_VDJ_CDR3CDR2_compare = shared_IgE_main.copy()
            
            # Update the comparison IgE DataFrame with boolean masks
            divergent_tab_low_moderate_high_mutated_shared_IgE_VDJ_CDR3CDR2_compare['Match_in_low_mutated_shared_IgE_VDJ_CDR3CDR2_2_seqid'] = divergent_tab_low_moderate_high_mutated_shared_IgE_VDJ_CDR3CDR2_compare['seqid'].isin(divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2_2_seqid)
            divergent_tab_low_moderate_high_mutated_shared_IgE_VDJ_CDR3CDR2_compare['Match_in_moderate_mutated_shared_IgE_VDJ_CDR3CDR2_2_seqid'] = divergent_tab_low_moderate_high_mutated_shared_IgE_VDJ_CDR3CDR2_compare['seqid'].isin(divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2_2_seqid)
            divergent_tab_low_moderate_high_mutated_shared_IgE_VDJ_CDR3CDR2_compare['Match_in_high_mutated_shared_IgE_VDJ_CDR3CDR2_2_seqid'] = divergent_tab_low_moderate_high_mutated_shared_IgE_VDJ_CDR3CDR2_compare['seqid'].isin(divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2_2_seqid)
            
            
            # Ensure 'Match_in_low_moderate and high mutated_shared_IgE_seqid' is a string and remove leading/trailing spaces
            divergent_tab_low_moderate_high_mutated_shared_IgE_VDJ_CDR3CDR2_compare['Match_in_low_mutated_shared_IgE_VDJ_CDR3CDR2_2_seqid'] = divergent_tab_low_moderate_high_mutated_shared_IgE_VDJ_CDR3CDR2_compare['Match_in_low_mutated_shared_IgE_VDJ_CDR3CDR2_2_seqid'].astype(str).str.strip()
            divergent_tab_low_moderate_high_mutated_shared_IgE_VDJ_CDR3CDR2_compare['Match_in_moderate_mutated_shared_IgE_VDJ_CDR3CDR2_2_seqid'] = divergent_tab_low_moderate_high_mutated_shared_IgE_VDJ_CDR3CDR2_compare['Match_in_moderate_mutated_shared_IgE_VDJ_CDR3CDR2_2_seqid'].astype(str).str.strip()
            divergent_tab_low_moderate_high_mutated_shared_IgE_VDJ_CDR3CDR2_compare['Match_in_high_mutated_shared_IgE_VDJ_CDR3CDR2_2_seqid'] = divergent_tab_low_moderate_high_mutated_shared_IgE_VDJ_CDR3CDR2_compare['Match_in_high_mutated_shared_IgE_VDJ_CDR3CDR2_2_seqid'].astype(str).str.strip()
            
            
            # Now filter the low_moderate and high mutated_shared_IgE clones
            divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone = divergent_tab_low_moderate_high_mutated_shared_IgE_VDJ_CDR3CDR2_compare[(divergent_tab_low_moderate_high_mutated_shared_IgE_VDJ_CDR3CDR2_compare['Match_in_low_mutated_shared_IgE_VDJ_CDR3CDR2_2_seqid'] == 'True')]
            divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone = divergent_tab_low_moderate_high_mutated_shared_IgE_VDJ_CDR3CDR2_compare[(divergent_tab_low_moderate_high_mutated_shared_IgE_VDJ_CDR3CDR2_compare['Match_in_moderate_mutated_shared_IgE_VDJ_CDR3CDR2_2_seqid'] == 'True')]
            divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone = divergent_tab_low_moderate_high_mutated_shared_IgE_VDJ_CDR3CDR2_compare[(divergent_tab_low_moderate_high_mutated_shared_IgE_VDJ_CDR3CDR2_compare['Match_in_high_mutated_shared_IgE_VDJ_CDR3CDR2_2_seqid'] == 'True')]
            
            # Determine the number of divergent low, moderate, and high mutated IgE clones
            unique_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone = divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone['IgE_VDJ'].nunique(dropna=True)
            print(f'Number of divergent CDR3CDR2 low mutated IgE clones: {unique_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone}')
            
            unique_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone = divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone['IgE_VDJ'].nunique(dropna=True)
            print(f'Number of divergent CDR3CDR2 moderate mutated IgE clones: {unique_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone}')
            
            unique_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone = divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone['IgE_VDJ'].nunique(dropna=True)
            print(f'Number of divergent CDR3CDR2 high mutated IgE clones: {unique_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone}')
            
            # Determine the average copy number of divergent low, moderate and high mutated IgE clones
            
            # Calculate the mean copynumber for low mutated IgE clones, ensuring NaN results are handled
            mean_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone =divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone['onecopy'].mean()
            mean_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone = 0 if pd.isna(mean_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone) else mean_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone
            print(f'mean copynumber of divergent CDR3CDR2 low mutated IgE: {mean_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone}')
            
            # Calculate the mean copynumber for moderate mutated IgE clones, ensuring NaN results are handled
            mean_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone = divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone['onecopy'].mean()
            mean_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone = 0 if pd.isna(mean_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone) else mean_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone
            print(f'mean copynumber of divergent CDR3CDR2 moderate mutated IgE: {mean_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone}')
            
            # Calculate the mean copynumber for high mutated IgE clones, ensuring NaN results are handled
            mean_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone = divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone['onecopy'].mean()
            mean_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone = 0 if pd.isna(mean_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone) else mean_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone
            print(f'mean copynumber of divergent CDR3CDR2 high mutated IgE: {mean_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone}')
            
            
            # Calculate the sum copynumber for low mutated IgE clones, ensuring NaN results are handled
            sum_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone =divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone['onecopy'].sum()
            sum_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone = 0 if pd.isna(sum_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone) else sum_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone
            print(f'sum copynumber of divergent CDR3CDR2 low mutated IgE: {sum_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone}')
            
            # Calculate the mean copynumber for moderate mutated IgE clones, ensuring NaN results are handled
            sum_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone = divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone['onecopy'].sum()
            sum_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone = 0 if pd.isna(sum_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone) else sum_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone
            print(f'sum copynumber of divergent CDR3CDR2 moderate mutated IgE: {sum_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone}')
            
            # Calculate the mean copynumber for high mutated IgE clones, ensuring NaN results are handled
            sum_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone = divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone['onecopy'].sum()
            sum_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone = 0 if pd.isna(sum_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone) else sum_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone
            print(f'sum copynumber of divergent CDR3CDR2 high mutated IgE: {sum_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone}')
            
            
            
            # Filter out the clonally related divergent low, moderate and high mutated CDR3CDR2 IgG1 counterpart
            # Create sets of low, moderate and high IgE CDR3.2.seqid for faster lookup
            
            divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone_seqid = set(divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone['IgE_VDJ'])
            divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone_seqid = set(divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone['IgE_VDJ'])
            divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone_seqid = set(divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone['IgE_VDJ'])
            
            # Create a copy of the data frame
            clonally_related_low_moderate_high_shared_CDR3CDR2_IgG1 = shared_IgG1_main.copy()
            
            # Update the comparison IgE DataFrame with boolean masks
            clonally_related_low_moderate_high_shared_CDR3CDR2_IgG1['divergent_low_mutated_shared_IgE_VDJ_CDR3CDR2_in_IgG1_seqid'] = clonally_related_low_moderate_high_shared_CDR3CDR2_IgG1['IgG1_VDJ'].isin(divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone_seqid)
            clonally_related_low_moderate_high_shared_CDR3CDR2_IgG1['divergent_moderate_mutated_shared_IgE_VDJ_CDR3CDR2_in_IgG1_seqid'] = clonally_related_low_moderate_high_shared_CDR3CDR2_IgG1['IgG1_VDJ'].isin(divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone_seqid)
            clonally_related_low_moderate_high_shared_CDR3CDR2_IgG1['divergent_high_mutated_shared_IgE_VDJ_CDR3CDR2_in_IgG1_seqid'] = clonally_related_low_moderate_high_shared_CDR3CDR2_IgG1['IgG1_VDJ'].isin(divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone_seqid)
            
            # Ensure 'Match_in_low_moderate and high mutated_shared_IgE_seqid' is a string and remove leading/trailing spaces
            clonally_related_low_moderate_high_shared_CDR3CDR2_IgG1['divergent_low_mutated_shared_IgE_VDJ_CDR3CDR2_in_IgG1_seqid'] = clonally_related_low_moderate_high_shared_CDR3CDR2_IgG1['divergent_low_mutated_shared_IgE_VDJ_CDR3CDR2_in_IgG1_seqid'].astype(str).str.strip()
            clonally_related_low_moderate_high_shared_CDR3CDR2_IgG1['divergent_moderate_mutated_shared_IgE_VDJ_CDR3CDR2_in_IgG1_seqid'] = clonally_related_low_moderate_high_shared_CDR3CDR2_IgG1['divergent_moderate_mutated_shared_IgE_VDJ_CDR3CDR2_in_IgG1_seqid'].astype(str).str.strip()
            clonally_related_low_moderate_high_shared_CDR3CDR2_IgG1['divergent_high_mutated_shared_IgE_VDJ_CDR3CDR2_in_IgG1_seqid'] = clonally_related_low_moderate_high_shared_CDR3CDR2_IgG1['divergent_high_mutated_shared_IgE_VDJ_CDR3CDR2_in_IgG1_seqid'].astype(str).str.strip()
            
            # Now filter the low_moderate and high mutated_shared_IgE clones
            clonally_related_low_mutated_shared_CDR3CDR2_IgG1 = clonally_related_low_moderate_high_shared_CDR3CDR2_IgG1[(clonally_related_low_moderate_high_shared_CDR3CDR2_IgG1['divergent_low_mutated_shared_IgE_VDJ_CDR3CDR2_in_IgG1_seqid'] == 'True')]
            clonally_related_moderate_mutated_shared_CDR3CDR2_IgG1 = clonally_related_low_moderate_high_shared_CDR3CDR2_IgG1[(clonally_related_low_moderate_high_shared_CDR3CDR2_IgG1['divergent_moderate_mutated_shared_IgE_VDJ_CDR3CDR2_in_IgG1_seqid'] == 'True')]
            clonally_related_high_mutated_shared_CDR3CDR2_IgG1 = clonally_related_low_moderate_high_shared_CDR3CDR2_IgG1[(clonally_related_low_moderate_high_shared_CDR3CDR2_IgG1['divergent_high_mutated_shared_IgE_VDJ_CDR3CDR2_in_IgG1_seqid'] == 'True')]
            
            # Determine the average copy number of divergent low, moderate and high mutated clonally related IgG1 counterpart
            # Calculate the mean copynumber for low mutated IgE clones IgG1 counterpart, ensuring NaN results are handled
            mean_clonally_related_low_mutated_shared_CDR3CDR2_IgG1 = clonally_related_low_mutated_shared_CDR3CDR2_IgG1['onecopy'].mean()
            # Handle NaN values explicitly
            mean_clonally_related_low_mutated_shared_CDR3CDR2_IgG1 = (
            0 if pd.isna(mean_clonally_related_low_mutated_shared_CDR3CDR2_IgG1)
            else mean_clonally_related_low_mutated_shared_CDR3CDR2_IgG1)
            print(f'mean copynumber of IgG1 clonally  related to low mutated IgE (CDR3CDR2): {mean_clonally_related_low_mutated_shared_CDR3CDR2_IgG1}')
            
            
            # Calculate the mean copynumber for moderate mutated IgE clones IgG1 counterpart, ensuring NaN results are handled
            mean_clonally_related_moderate_mutated_shared_CDR3CDR2_IgG1 = clonally_related_moderate_mutated_shared_CDR3CDR2_IgG1['onecopy'].mean()
            mean_clonally_related_moderate_mutated_shared_CDR3CDR2_IgG1 = 0 if pd.isna(mean_clonally_related_moderate_mutated_shared_CDR3CDR2_IgG1) else mean_clonally_related_moderate_mutated_shared_CDR3CDR2_IgG1
            print(f'mean copynumber of IgG1 clonally related to moderate mutated IgE (CDR3CDR2): {mean_clonally_related_moderate_mutated_shared_CDR3CDR2_IgG1}')
            
            # Calculate the mean copynumber for moderate mutated IgE clones IgG1 counterpart, ensuring NaN results are handled
            mean_clonally_related_high_mutated_shared_CDR3CDR2_IgG1 = clonally_related_high_mutated_shared_CDR3CDR2_IgG1['onecopy'].mean()
            mean_clonally_related_high_mutated_shared_CDR3CDR2_IgG1 = 0 if pd.isna(mean_clonally_related_high_mutated_shared_CDR3CDR2_IgG1) else mean_clonally_related_high_mutated_shared_CDR3CDR2_IgG1
            print(f'mean copynumber of IgG1 clonally related to high mutated IgE (CDR3CDR2): {mean_clonally_related_high_mutated_shared_CDR3CDR2_IgG1}')
            
            # Filter out the entire VH gene of the low, moderate and highly mutated IgE using CDR3 AA change to filter the divergent IgEs
            #make a copy of the dataframe
            stats_IgE_Aminotab_productive_CDR3mut2 = stats_IgE_Aminotab_productive.copy()
            
            # Update the comparison IgE DataFrame with boolean masks
            stats_IgE_Aminotab_productive_CDR3mut2['divergent_low_mutated_shared_IgE_VH_AA'] = stats_IgE_Aminotab_productive_CDR3mut2['Sequence ID'].isin(divergent_tab_low_mutated_shared_IgE_VDJ_CDR3_2_seqid)
            stats_IgE_Aminotab_productive_CDR3mut2['divergent_moderate_mutated_shared_IgE_VH_AA'] = stats_IgE_Aminotab_productive_CDR3mut2['Sequence ID'].isin(divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3_2_seqid)
            stats_IgE_Aminotab_productive_CDR3mut2['divergent_high_mutated_shared_IgE_VH_AA'] = stats_IgE_Aminotab_productive_CDR3mut2['Sequence ID'].isin(divergent_tab_high_mutated_shared_IgE_VDJ_CDR3_2_seqid)
            
            # Ensure 'Match_in_low_moderate and high mutated_shared_IgE_seqid' is a string and remove leading/trailing spaces
            stats_IgE_Aminotab_productive_CDR3mut2['divergent_low_mutated_shared_IgE_VH_AA'] = stats_IgE_Aminotab_productive_CDR3mut2['divergent_low_mutated_shared_IgE_VH_AA'].astype(str).str.strip()
            stats_IgE_Aminotab_productive_CDR3mut2['divergent_moderate_mutated_shared_IgE_VH_AA'] = stats_IgE_Aminotab_productive_CDR3mut2['divergent_moderate_mutated_shared_IgE_VH_AA'].astype(str).str.strip()
            stats_IgE_Aminotab_productive_CDR3mut2['divergent_high_mutated_shared_IgE_VH_AA'] = stats_IgE_Aminotab_productive_CDR3mut2['divergent_high_mutated_shared_IgE_VH_AA'].astype(str).str.strip()
            
            # Now filter the low_moderate and high mutated_shared_IgE clones
            divergent_low_mutated_shared_IgE_VH_AA_tab = stats_IgE_Aminotab_productive_CDR3mut2[(stats_IgE_Aminotab_productive_CDR3mut2['divergent_low_mutated_shared_IgE_VH_AA'] == 'True')]
            divergent_moderate_mutated_shared_IgE_VH_AA_tab = stats_IgE_Aminotab_productive_CDR3mut2[(stats_IgE_Aminotab_productive_CDR3mut2['divergent_moderate_mutated_shared_IgE_VH_AA'] == 'True')]
            divergent_high_mutated_shared_IgE_VH_AA_tab = stats_IgE_Aminotab_productive_CDR3mut2[(stats_IgE_Aminotab_productive_CDR3mut2['divergent_high_mutated_shared_IgE_VH_AA'] == 'True')]
            
            
            # Extract sequences from V-REGION column for low mutated V region
            divergent_low_mutated_VH_AAsequences = divergent_low_mutated_shared_IgE_VH_AA_tab["V-D-J-REGION"].tolist()
            
            # Print the extracted sequences
            for i, sequence in enumerate(divergent_low_mutated_VH_AAsequences, 1):
                print(f"Row {i}: {sequence}")
            
            
            # Extract sequences from V-REGION column for moderate mutated V region
            divergent_moderate_mutated_VH_AAsequences = divergent_moderate_mutated_shared_IgE_VH_AA_tab["V-D-J-REGION"].tolist()
            
            # Print the extracted sequences
            for i, sequence in enumerate(divergent_moderate_mutated_VH_AAsequences, 1):
                print(f"Row {i}: {sequence}")
            
            
            # Extract sequences from V-REGION column for high mutated V region
            divergent_high_mutated_VH_AAsequences = divergent_high_mutated_shared_IgE_VH_AA_tab["V-D-J-REGION"].tolist()
            
            # Print the extracted sequences
            for i, sequence in enumerate(divergent_high_mutated_VH_AAsequences, 1):
                print(f"Row {i}: {sequence}")
            
            
            
            # Filter out the entire VH gene of the clonally related IgG1 to the low, moderate and highly mutated IgE
            clonally_related_low_mutated_shared_IgG1_seqid = set(clonally_related_low_mutated_shared_IgG1['seqid'])
            clonally_related_moderate_mutated_shared_IgG1_seqid = set(clonally_related_moderate_mutated_shared_IgG1['seqid'])
            clonally_related_high_mutated_shared_IgG1_seqid = set(clonally_related_high_mutated_shared_IgG1['seqid'])
            
            # Create a copy of the data frame
            stats_IgG1_Aminotab_productive_CDR3mut2 = stats_IgG1_Aminotab_productive.copy()
            
            # Update the comparison IgG1 DataFrame with boolean masks
            # Check for the correct column name in stats_IgG1_Aminotab_productive_CDR3mut2
            if 'Sequence.ID' in stats_IgG1_Aminotab_productive_CDR3mut2.columns:
                column_seqid = 'Sequence.ID'
            elif 'Sequence ID' in stats_IgG1_Aminotab_productive_CDR3mut2.columns:
                column_seqid = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present in the DataFrame.")
            
            # Perform the analysis by checking if sequence IDs are in clonally_related_low_mutated_shared_IgG1_seqid
            stats_IgG1_Aminotab_productive_CDR3mut2['clonally_related_low_mutated_shared_IgG1_seqid_in_IgG1_seqid'] = \
                stats_IgG1_Aminotab_productive_CDR3mut2[column_seqid].isin(clonally_related_low_mutated_shared_IgG1_seqid)


            # Check for the correct column name in stats_IgG1_Aminotab_productive_CDR3mut2
            if 'Sequence.ID' in stats_IgG1_Aminotab_productive_CDR3mut2.columns:
                column_seqid = 'Sequence.ID'
            elif 'Sequence ID' in stats_IgG1_Aminotab_productive_CDR3mut2.columns:
                column_seqid = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present in the DataFrame.")
            
            # Perform the analysis by checking if sequence IDs are in clonally_related_moderate_mutated_shared_IgG1_seqid
            stats_IgG1_Aminotab_productive_CDR3mut2['clonally_related_moderate_mutated_shared_IgG1_seqid_in_IgG1_seqid'] = \
                stats_IgG1_Aminotab_productive_CDR3mut2[column_seqid].isin(clonally_related_moderate_mutated_shared_IgG1_seqid)


            # Check for the correct column name in stats_IgG1_Aminotab_productive_CDR3mut2
            if 'Sequence.ID' in stats_IgG1_Aminotab_productive_CDR3mut2.columns:
                column_seqid = 'Sequence.ID'
            elif 'Sequence ID' in stats_IgG1_Aminotab_productive_CDR3mut2.columns:
                column_seqid = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present in the DataFrame.")
            
            # Perform the analysis by checking if sequence IDs are in clonally_related_high_mutated_shared_IgG1_seqid
            stats_IgG1_Aminotab_productive_CDR3mut2['clonally_related_high_mutated_shared_IgG1_seqid_in_IgG1_seqid'] = \
                stats_IgG1_Aminotab_productive_CDR3mut2[column_seqid].isin(clonally_related_high_mutated_shared_IgG1_seqid)

            
            # Ensure 'Match_in_low_moderate and high mutated_shared_IgE_seqid' is a string and remove leading/trailing spaces
            stats_IgG1_Aminotab_productive_CDR3mut2['clonally_related_low_mutated_shared_IgG1_seqid_in_IgG1_seqid'] = stats_IgG1_Aminotab_productive_CDR3mut2['clonally_related_low_mutated_shared_IgG1_seqid_in_IgG1_seqid'].astype(str).str.strip()
            stats_IgG1_Aminotab_productive_CDR3mut2['clonally_related_moderate_mutated_shared_IgG1_seqid_in_IgG1_seqid'] = stats_IgG1_Aminotab_productive_CDR3mut2['clonally_related_moderate_mutated_shared_IgG1_seqid_in_IgG1_seqid'].astype(str).str.strip()
            stats_IgG1_Aminotab_productive_CDR3mut2['clonally_related_high_mutated_shared_IgG1_seqid_in_IgG1_seqid'] = stats_IgG1_Aminotab_productive_CDR3mut2['clonally_related_high_mutated_shared_IgG1_seqid_in_IgG1_seqid'].astype(str).str.strip()
            
            # Now filter the low_moderate and high mutated_shared_clonally related IgG1 clones
            clonally_related_low_mutated_shared_IgG1_seqid_in_IgG1_AA = stats_IgG1_Aminotab_productive_CDR3mut2[(stats_IgG1_Aminotab_productive_CDR3mut2['clonally_related_low_mutated_shared_IgG1_seqid_in_IgG1_seqid'] == 'True')]
            clonally_related_moderate_mutated_shared_IgG1_seqid_in_IgG1_AA = stats_IgG1_Aminotab_productive_CDR3mut2[(stats_IgG1_Aminotab_productive_CDR3mut2['clonally_related_moderate_mutated_shared_IgG1_seqid_in_IgG1_seqid'] == 'True')]
            clonally_related_high_mutated_shared_IgG1_seqid_in_IgG1_AA = stats_IgG1_Aminotab_productive_CDR3mut2[(stats_IgG1_Aminotab_productive_CDR3mut2['clonally_related_high_mutated_shared_IgG1_seqid_in_IgG1_seqid'] == 'True')]
            
            
            # Extract sequences from V-REGION column for low mutated V region

            # Check for the correct column name in clonally_related_low_mutated_shared_IgG1_seqid_in_IgG1_AA
            if 'V.D.J.REGION' in clonally_related_low_mutated_shared_IgG1_seqid_in_IgG1_AA.columns:
                column_vdj_region = 'V.D.J.REGION'
            elif 'V-D-J-REGION' in clonally_related_low_mutated_shared_IgG1_seqid_in_IgG1_AA.columns:
                column_vdj_region = 'V-D-J-REGION'
            else:
                raise KeyError("Neither 'V.D.J.REGION' nor 'V-D-J-REGION' column is present in the DataFrame.")
            
            # Convert the selected column to a list
            clonally_related_low_mutated_shared_IgG1_seqid_in_IgG1_AAsequences = \
                clonally_related_low_mutated_shared_IgG1_seqid_in_IgG1_AA[column_vdj_region].tolist()

            
            # Print the extracted sequences
            for i, sequence in enumerate(clonally_related_low_mutated_shared_IgG1_seqid_in_IgG1_AAsequences, 1):
                print(f"Row {i}: {sequence}")
            
            
            # Extract sequences from V-REGION column for moderate mutated V region
            
            # Check for the correct column name in clonally_related_moderate_mutated_shared_IgG1_seqid_in_IgG1_AA
            if 'V.D.J.REGION' in clonally_related_moderate_mutated_shared_IgG1_seqid_in_IgG1_AA.columns:
                column_vdj_region = 'V.D.J.REGION'
            elif 'V-D-J-REGION' in clonally_related_moderate_mutated_shared_IgG1_seqid_in_IgG1_AA.columns:
                column_vdj_region = 'V-D-J-REGION'
            else:
                raise KeyError("Neither 'V.D.J.REGION' nor 'V-D-J-REGION' column is present in the DataFrame.")
            
            # Convert the selected column to a list
            clonally_related_moderate_mutated_shared_IgG1_seqid_in_IgG1_AAsequences = \
                clonally_related_moderate_mutated_shared_IgG1_seqid_in_IgG1_AA[column_vdj_region].tolist()

            
            # Print the extracted sequences
            for i, sequence in enumerate(clonally_related_moderate_mutated_shared_IgG1_seqid_in_IgG1_AAsequences, 1):
                print(f"Row {i}: {sequence}")
            
            
            # Extract sequences from V-REGION column for high mutated V region

            # Check for the correct column name in clonally_related_high_mutated_shared_IgG1_seqid_in_IgG1_AA
            if 'V.D.J.REGION' in clonally_related_high_mutated_shared_IgG1_seqid_in_IgG1_AA.columns:
                column_vdj_region = 'V.D.J.REGION'
            elif 'V-D-J-REGION' in clonally_related_high_mutated_shared_IgG1_seqid_in_IgG1_AA.columns:
                column_vdj_region = 'V-D-J-REGION'
            else:
                raise KeyError("Neither 'V.D.J.REGION' nor 'V-D-J-REGION' column is present in the DataFrame.")
            
            # Convert the selected column to a list
            clonally_related_high_mutated_shared_IgG1_seqid_in_IgG1_AAsequences = \
                clonally_related_high_mutated_shared_IgG1_seqid_in_IgG1_AA[column_vdj_region].tolist()


            
            # Print the extracted sequences
            for i, sequence in enumerate(clonally_related_high_mutated_shared_IgG1_seqid_in_IgG1_AAsequences, 1):
                print(f"Row {i}: {sequence}")
            
            
            
            # Filter out the entire VH gene of the low, moderate and highly mutated IgE using CDR3 AA change to filter the non-divergent IgEs
            #make a copy of the dataframe
            stats_IgE_Aminotab_productive_CDR3mut2 = stats_IgE_Aminotab_productive.copy()
            
            # Update the comparison IgE DataFrame with boolean masks
            stats_IgE_Aminotab_productive_CDR3mut2['non_divergent_low_mutated_shared_IgE_VH_AA'] = stats_IgE_Aminotab_productive_CDR3mut2['Sequence ID'].isin(non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3_2_seqid)
            stats_IgE_Aminotab_productive_CDR3mut2['non_divergent_moderate_mutated_shared_IgE_VH_AA'] = stats_IgE_Aminotab_productive_CDR3mut2['Sequence ID'].isin(non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3_2_seqid)
            stats_IgE_Aminotab_productive_CDR3mut2['non_divergent_high_mutated_shared_IgE_VH_AA'] = stats_IgE_Aminotab_productive_CDR3mut2['Sequence ID'].isin(non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3_2_seqid)
            
            # Ensure 'Match_in_low_moderate and high mutated_shared_IgE_seqid' is a string and remove leading/trailing spaces
            stats_IgE_Aminotab_productive_CDR3mut2['non_divergent_low_mutated_shared_IgE_VH_AA'] = stats_IgE_Aminotab_productive_CDR3mut2['non_divergent_low_mutated_shared_IgE_VH_AA'].astype(str).str.strip()
            stats_IgE_Aminotab_productive_CDR3mut2['non_divergent_moderate_mutated_shared_IgE_VH_AA'] = stats_IgE_Aminotab_productive_CDR3mut2['non_divergent_moderate_mutated_shared_IgE_VH_AA'].astype(str).str.strip()
            stats_IgE_Aminotab_productive_CDR3mut2['non_divergent_high_mutated_shared_IgE_VH_AA'] = stats_IgE_Aminotab_productive_CDR3mut2['non_divergent_high_mutated_shared_IgE_VH_AA'].astype(str).str.strip()
            
            # Now filter the low_moderate and high mutated_shared_IgE clones
            non_divergent_low_mutated_shared_IgE_VH_AA_tab = stats_IgE_Aminotab_productive_CDR3mut2[(stats_IgE_Aminotab_productive_CDR3mut2['non_divergent_low_mutated_shared_IgE_VH_AA'] == 'True')]
            non_divergent_moderate_mutated_shared_IgE_VH_AA_tab = stats_IgE_Aminotab_productive_CDR3mut2[(stats_IgE_Aminotab_productive_CDR3mut2['non_divergent_moderate_mutated_shared_IgE_VH_AA'] == 'True')]
            non_divergent_high_mutated_shared_IgE_VH_AA_tab = stats_IgE_Aminotab_productive_CDR3mut2[(stats_IgE_Aminotab_productive_CDR3mut2['non_divergent_high_mutated_shared_IgE_VH_AA'] == 'True')]
            
            
            # Extract sequences from V-REGION column for low mutated V region
            non_divergent_low_mutated_VH_AAsequences = non_divergent_low_mutated_shared_IgE_VH_AA_tab["V-D-J-REGION"].tolist()
            
            # Print the extracted sequences
            for i, sequence in enumerate(non_divergent_low_mutated_VH_AAsequences, 1):
                print(f"Row {i}: {sequence}")
            
            
            # Extract sequences from V-REGION column for moderate mutated V region
            non_divergent_moderate_mutated_VH_AAsequences = non_divergent_moderate_mutated_shared_IgE_VH_AA_tab["V-D-J-REGION"].tolist()
            
            # Print the extracted sequences
            for i, sequence in enumerate(non_divergent_moderate_mutated_VH_AAsequences, 1):
                print(f"Row {i}: {sequence}")
            
            
            # Extract sequences from V-REGION column for high mutated V region
            non_divergent_high_mutated_VH_AAsequences = non_divergent_high_mutated_shared_IgE_VH_AA_tab["V-D-J-REGION"].tolist()
            
            # Print the extracted sequences
            for i, sequence in enumerate(non_divergent_high_mutated_VH_AAsequences, 1):
                print(f"Row {i}: {sequence}")
            
            
            # Filter out the entire VH gene of the top ten low, moderate and highly mutated IgE
            # Make a copy of the original DataFrame
            stats_IgE_Aminotab_productive_CDR3mut2 = stats_IgE_Aminotab_productive.copy()
            
            # Create sets of top ten low, moderate, and high mutated shared IgE seqids
            stats_tab_low_mutated_shared_IgE_VDJ_sorted_topten_seqid = set(stats_tab_low_mutated_shared_IgE_VDJ_sorted_topten['seqid'])
            stats_tab_moderate_mutated_shared_IgE_VDJ_sorted_topten_seqid = set(stats_tab_moderate_mutated_shared_IgE_VDJ_sorted_topten['seqid'])
            stats_tab_high_mutated_shared_IgE_VDJ_sorted_topten_seqid = set(stats_tab_high_mutated_shared_IgE_VDJ_sorted_topten['seqid'])
            
            # Add boolean columns for membership in top-ten sets
            stats_IgE_Aminotab_productive_CDR3mut2['topten_low_mutated_shared_IgE_VH_AA'] = stats_IgE_Aminotab_productive_CDR3mut2['Sequence ID'].isin(stats_tab_low_mutated_shared_IgE_VDJ_sorted_topten_seqid)
            stats_IgE_Aminotab_productive_CDR3mut2['topten_moderate_mutated_shared_IgE_VH_AA'] = stats_IgE_Aminotab_productive_CDR3mut2['Sequence ID'].isin(stats_tab_moderate_mutated_shared_IgE_VDJ_sorted_topten_seqid)
            stats_IgE_Aminotab_productive_CDR3mut2['topten_high_mutated_shared_IgE_VH_AA'] = stats_IgE_Aminotab_productive_CDR3mut2['Sequence ID'].isin(stats_tab_high_mutated_shared_IgE_VDJ_sorted_topten_seqid)
            
            # Filter rows based on top-ten membership
            stats_tab_low_mutated_shared_IgE_VDJ_sorted_topten_tab = stats_IgE_Aminotab_productive_CDR3mut2[stats_IgE_Aminotab_productive_CDR3mut2['topten_low_mutated_shared_IgE_VH_AA']]
            stats_tab_moderate_mutated_shared_IgE_VDJ_sorted_topten_tab = stats_IgE_Aminotab_productive_CDR3mut2[stats_IgE_Aminotab_productive_CDR3mut2['topten_moderate_mutated_shared_IgE_VH_AA']]
            stats_tab_high_mutated_shared_IgE_VDJ_sorted_topten_tab = stats_IgE_Aminotab_productive_CDR3mut2[stats_IgE_Aminotab_productive_CDR3mut2['topten_high_mutated_shared_IgE_VH_AA']]
            
            # Extract sequences from the V-D-J-REGION column for each category
            stats_tab_low_mutated_shared_IgE_VDJ_sorted_topten_AA = stats_tab_low_mutated_shared_IgE_VDJ_sorted_topten_tab["V-D-J-REGION"].tolist()
            stats_tab_moderate_mutated_shared_IgE_VDJ_sorted_topten_AA = stats_tab_moderate_mutated_shared_IgE_VDJ_sorted_topten_tab["V-D-J-REGION"].tolist()
            stats_tab_high_mutated_shared_IgE_VDJ_sorted_topten_AA = stats_tab_high_mutated_shared_IgE_VDJ_sorted_topten_tab["V-D-J-REGION"].tolist()
            
            # Print extracted sequences for each category
            print("\n top ten Low Mutated Shared IgE Sequences:")
            for i, sequence in enumerate(stats_tab_low_mutated_shared_IgE_VDJ_sorted_topten_AA, 1):
                print(f"Row {i} (Length {len(sequence)}): {sequence}")
            
            print("\n top ten Moderate Mutated Shared IgE Sequences:")
            for i, sequence in enumerate(stats_tab_moderate_mutated_shared_IgE_VDJ_sorted_topten_AA, 1):
                print(f"Row {i} (Length {len(sequence)}): {sequence}")
            
            print("\n top ten High Mutated Shared IgE Sequences:")
            for i, sequence in enumerate(stats_tab_high_mutated_shared_IgE_VDJ_sorted_topten_AA, 1):
                print(f"Row {i} (Length {len(sequence)}): {sequence}")
            
            
            
            # Filter out the entire VH gene of the top ten low, moderate and highly mutated IgG1
            
            # Make a copy of the original DataFrame
            stats_IgG1_Aminotab_productive_CDR3mut2 = stats_IgG1_Aminotab_productive.copy()
            
            # Create sets of top ten low, moderate, and high mutated shared IgE seqids
            stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten_seqid = set(stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten['seqid'])
            stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten_seqid = set(stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten['seqid'])
            stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten_seqid = set(stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten['seqid'])
            
            # Add boolean columns for membership in top-ten sets
            # Check for the correct column name in stats_IgG1_Aminotab_productive_CDR3mut2
            if 'Sequence.ID' in stats_IgG1_Aminotab_productive_CDR3mut2.columns:
                column_seqid = 'Sequence.ID'
            elif 'Sequence ID' in stats_IgG1_Aminotab_productive_CDR3mut2.columns:
                column_seqid = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present in the DataFrame.")
            
            # Perform the analysis by checking if sequence IDs are in stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten_seqid
            stats_IgG1_Aminotab_productive_CDR3mut2['topten_low_mutated_shared_IgG1_VH_AA'] = \
                stats_IgG1_Aminotab_productive_CDR3mut2[column_seqid].isin(stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten_seqid)


            # Check for the correct column name in stats_IgG1_Aminotab_productive_CDR3mut2
            if 'Sequence.ID' in stats_IgG1_Aminotab_productive_CDR3mut2.columns:
                column_seqid = 'Sequence.ID'
            elif 'Sequence ID' in stats_IgG1_Aminotab_productive_CDR3mut2.columns:
                column_seqid = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present in the DataFrame.")
            
            # Perform the analysis by checking if sequence IDs are in stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten_seqid
            stats_IgG1_Aminotab_productive_CDR3mut2['topten_moderate_mutated_shared_IgG1_VH_AA'] = \
                stats_IgG1_Aminotab_productive_CDR3mut2[column_seqid].isin(stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten_seqid)



            # Check for the correct column name in stats_IgG1_Aminotab_productive_CDR3mut2
            if 'Sequence.ID' in stats_IgG1_Aminotab_productive_CDR3mut2.columns:
                column_seqid = 'Sequence.ID'
            elif 'Sequence ID' in stats_IgG1_Aminotab_productive_CDR3mut2.columns:
                column_seqid = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present in the DataFrame.")
            
            # Perform the analysis by checking if sequence IDs are in stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten_seqid
            stats_IgG1_Aminotab_productive_CDR3mut2['topten_high_mutated_shared_IgG1_VH_AA'] = \
                stats_IgG1_Aminotab_productive_CDR3mut2[column_seqid].isin(stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten_seqid)

            
            # Filter rows based on top-ten membership
            stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten_tab = stats_IgG1_Aminotab_productive_CDR3mut2[stats_IgG1_Aminotab_productive_CDR3mut2['topten_low_mutated_shared_IgG1_VH_AA']]
            stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten_tab = stats_IgG1_Aminotab_productive_CDR3mut2[stats_IgG1_Aminotab_productive_CDR3mut2['topten_moderate_mutated_shared_IgG1_VH_AA']]
            stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten_tab = stats_IgG1_Aminotab_productive_CDR3mut2[stats_IgG1_Aminotab_productive_CDR3mut2['topten_high_mutated_shared_IgG1_VH_AA']]
            
            # Extract sequences from the V-D-J-REGION column for each category

            # Check for the correct column name in stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten_tab
            if 'V.D.J.REGION' in stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten_tab.columns:
                column_vdj_region = 'V.D.J.REGION'
            elif 'V-D-J-REGION' in stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten_tab.columns:
                column_vdj_region = 'V-D-J-REGION'
            else:
                raise KeyError("Neither 'V.D.J.REGION' nor 'V-D-J-REGION' column is present in the DataFrame.")
            
            # Convert the selected column to a list
            stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten_AA = \
                stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten_tab[column_vdj_region].tolist()


            # Check for the correct column name in stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten_tab
            if 'V.D.J.REGION' in stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten_tab.columns:
                column_vdj_region = 'V.D.J.REGION'
            elif 'V-D-J-REGION' in stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten_tab.columns:
                column_vdj_region = 'V-D-J-REGION'
            else:
                raise KeyError("Neither 'V.D.J.REGION' nor 'V-D-J-REGION' column is present in the DataFrame.")
            
            # Convert the selected column to a list
            stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten_AA = \
                stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten_tab[column_vdj_region].tolist()


            # Check for the correct column name in stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten_tab
            if 'V.D.J.REGION' in stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten_tab.columns:
                column_vdj_region = 'V.D.J.REGION'
            elif 'V-D-J-REGION' in stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten_tab.columns:
                column_vdj_region = 'V-D-J-REGION'
            else:
                raise KeyError("Neither 'V.D.J.REGION' nor 'V-D-J-REGION' column is present in the DataFrame.")
            
            # Convert the selected column to a list
            stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten_AA = \
                stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten_tab[column_vdj_region].tolist()
            
            # Print extracted sequences for each category
            print("\n top ten Low Mutated Shared IgG1 Sequences:")
            for i, sequence in enumerate(stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten_AA, 1):
                print(f"Row {i} (Length {len(sequence)}): {sequence}")
            
            print("\n top ten Moderate Mutated Shared IgG1 Sequences:")
            for i, sequence in enumerate(stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten_AA, 1):
                print(f"Row {i} (Length {len(sequence)}): {sequence}")
            
            print("\n top ten High Mutated Shared IgG1 Sequences:")
            for i, sequence in enumerate(stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten_AA, 1):
                print(f"Row {i} (Length {len(sequence)}): {sequence}")
            
            
            #Process (low,moderate and high) divergent IgE clones and procees top ten low,moderate and high mutated Shared IgE clones
            import os
            from pathlib import Path
            from Bio import SeqIO
            from Bio.Align import PairwiseAligner
            import pandas as pd
            
            # Load FASTA sequences
            def load_imgt_fasta(fasta_file):
                sequences = {}
                for record in SeqIO.parse(fasta_file, "fasta"):
                    sequences[record.id] = str(record.seq)
                return sequences
            
            # Find the best alignment
            def find_best_alignment(fasta_sequences, query_sequence):
                aligner = PairwiseAligner()
                aligner.mode = "global"
                best_score = -float('inf')
                best_reference = None
                best_seq_id = None
            
                for seq_id, reference_sequence in fasta_sequences.items():
                    score = aligner.score(reference_sequence, query_sequence)
                    if score > best_score:
                        best_score = score
                        best_reference = reference_sequence
                        best_seq_id = seq_id
                
                return best_seq_id, best_reference
            
            # Calculate similarity between two sequences
            def calculate_similarity(seq1, seq2):
                matches = sum(a == b for a, b in zip(seq1, seq2))
                return matches / len(seq1)
            
            # Extract and concatenate sequences
            def extract_and_concatenate(query_sequence, reference_sequence):
                aligner = PairwiseAligner()
                aligner.mode = "global"
                alignment = aligner.align(reference_sequence, query_sequence)[0]
                aligned_reference = alignment.aligned[1]
            
                if len(aligned_reference) == 0:
                    print("No alignment found.")
                    return query_sequence
            
                start_index_reference = aligned_reference[-1][1]
                upstream_region = reference_sequence[:start_index_reference]
            
                overlap_length = 0
                for i in range(1, len(upstream_region) + 1):
                    overlap_candidate = upstream_region[-i:]
                    query_start_candidate = query_sequence[:i]
                    if len(overlap_candidate) == len(query_start_candidate):
                        similarity = calculate_similarity(overlap_candidate, query_start_candidate)
                        if similarity > 0.6:
                            overlap_length = i
                            break
                
                trimmed_upstream = upstream_region[:-overlap_length] if overlap_length > 0 else upstream_region
                return trimmed_upstream + query_sequence
            
            # Remove repeats in the sequence
            def keep_one_repeat(sequence, min_length=5):
                substrings = {}
                result = list(sequence)
            
                for i in range(len(sequence)):
                    for j in range(i + min_length, len(sequence) + 1):
                        substring = sequence[i:j]
                        if substring in substrings:
                            start_index = sequence.find(substring, substrings[substring] + 1)
                            if start_index != -1:
                                for k in range(start_index, start_index + len(substring)):
                                    result[k] = ""
                        else:
                            substrings[substring] = i
            
                return "".join(result)
            
            # Function to process and save completed sequences
            def process_and_save_sequences(input_df, fasta_sequences, output_filename):
                completed_sequences = []
                
                for _, row in input_df.iterrows():
                    query_sequence = row["V-D-J-REGION"]
                    _, best_reference = find_best_alignment(fasta_sequences, query_sequence)
            
                    completed_sequence = extract_and_concatenate(query_sequence, best_reference) if best_reference else query_sequence
                    completed_sequences.append(completed_sequence)
            
                # Add completed sequences and remove repeats
                input_df["Completed_VH_Sequences"] = completed_sequences
                input_df["cleaned_sequence"] = input_df["Completed_VH_Sequences"].apply(keep_one_repeat)
            
                # Get the path to the Downloads folder and save the output CSV
                downloads_folder = str(Path.home() / "Downloads")
                output_csv = os.path.join(downloads_folder, output_filename)
                input_df.to_csv(output_csv, index=False)
                print(f"Output saved to {output_csv}")
            
            # Main function to process different groups
            def main():
                fasta_file_path = self.file_paths[8]  # Replace with actual path
                fasta_sequences = load_imgt_fasta(fasta_file_path)
            
                # Process low mutated divergent shared IgE sequences
                process_and_save_sequences(
                    divergent_low_mutated_shared_IgE_VH_AA_tab.copy(),
                    fasta_sequences,
                    "output_divergent_low_mutated_shared_IgE_VH_AA_tab.csv"
                )
            
                # Process moderate mutated divergent shared IgE sequences
                process_and_save_sequences(
                    divergent_moderate_mutated_shared_IgE_VH_AA_tab.copy(),
                    fasta_sequences,
                    "output_divergent_moderate_mutated_shared_IgE_VH_AA_tab.csv"
                )
            
                # Process high mutated divergent shared IgE sequences
                process_and_save_sequences(
                    divergent_high_mutated_shared_IgE_VH_AA_tab.copy(),
                    fasta_sequences,
                    "output_divergent_high_mutated_shared_IgE_VH_AA_tab.csv"
                )
                
                # Process low mutated non-divergent shared IgE sequences
                process_and_save_sequences(
                    non_divergent_low_mutated_shared_IgE_VH_AA_tab.copy(),
                    fasta_sequences,
                    "output_non_divergent_low_mutated_shared_IgE_VH_AA_tab.csv"
                )
                
                # Process moderate mutated non-divergent shared IgE sequences
                process_and_save_sequences(
                    non_divergent_moderate_mutated_shared_IgE_VH_AA_tab.copy(),
                    fasta_sequences,
                    "output_non_divergent_moderate_mutated_shared_IgE_VH_AA_tab.csv"
                )
                
                # Process high mutated non-divergent shared IgE sequences
                process_and_save_sequences(
                    non_divergent_high_mutated_shared_IgE_VH_AA_tab.copy(),
                    fasta_sequences,
                    "output_non_divergent_high_mutated_shared_IgE_VH_AA_tab.csv"
                )
                
                # Process top ten low mutated shared IgE
                process_and_save_sequences(
                    stats_tab_low_mutated_shared_IgE_VDJ_sorted_topten_tab.copy(),
                    fasta_sequences,
                    "output_stats_tab_low_mutated_shared_IgE_VDJ_sorted_topten.csv" 
                )
                
               # Process top ten moderate mutated shared IgE
                process_and_save_sequences(
                    stats_tab_moderate_mutated_shared_IgE_VDJ_sorted_topten_tab.copy(),
                    fasta_sequences,
                    "output_stats_tab_moderate_mutated_shared_IgE_VDJ_sorted_topten.csv"
                )
                
                # Process top ten high mutated shared IgE
                process_and_save_sequences(
                    stats_tab_high_mutated_shared_IgE_VDJ_sorted_topten_tab.copy(),
                    fasta_sequences,
                    "output_stats_tab_high_mutated_shared_IgE_VDJ_sorted_topten.csv"
                )
            
            if __name__ == "__main__":
                main()
            
            #Analyse IgG1 clones clonally related to low, moderate and high divergent IgE, also process top ten low, moderate and high mutated shared IgG1 
            # Load FASTA sequences
            def load_imgt_fasta(fasta_file):
                sequences = {}
                for record in SeqIO.parse(fasta_file, "fasta"):
                    sequences[record.id] = str(record.seq)
                return sequences
            
            # Find the best alignment
            def find_best_alignment(fasta_sequences, query_sequence):
                aligner = PairwiseAligner()
                aligner.mode = "global"
                best_score = -float('inf')
                best_reference = None
                best_seq_id = None
            
                for seq_id, reference_sequence in fasta_sequences.items():
                    score = aligner.score(reference_sequence, query_sequence)
                    if score > best_score:
                        best_score = score
                        best_reference = reference_sequence
                        best_seq_id = seq_id
                
                return best_seq_id, best_reference
            
            # Calculate similarity between two sequences
            def calculate_similarity(seq1, seq2):
                matches = sum(a == b for a, b in zip(seq1, seq2))
                return matches / len(seq1)
            
            # Extract and concatenate sequences
            def extract_and_concatenate(query_sequence, reference_sequence):
                aligner = PairwiseAligner()
                aligner.mode = "global"
                alignment = aligner.align(reference_sequence, query_sequence)[0]
                aligned_reference = alignment.aligned[1]
            
                if len(aligned_reference) == 0:
                    print("No alignment found.")
                    return query_sequence
            
                start_index_reference = aligned_reference[-1][1]
                upstream_region = reference_sequence[:start_index_reference]
            
                overlap_length = 0
                for i in range(1, len(upstream_region) + 1):
                    overlap_candidate = upstream_region[-i:]
                    query_start_candidate = query_sequence[:i]
                    if len(overlap_candidate) == len(query_start_candidate):
                        similarity = calculate_similarity(overlap_candidate, query_start_candidate)
                        if similarity > 0.6:
                            overlap_length = i
                            break
                
                trimmed_upstream = upstream_region[:-overlap_length] if overlap_length > 0 else upstream_region
                return trimmed_upstream + query_sequence
            
            # Remove repeats in the sequence
            def keep_one_repeat(sequence, min_length=5):
                substrings = {}
                result = list(sequence)
            
                for i in range(len(sequence)):
                    for j in range(i + min_length, len(sequence) + 1):
                        substring = sequence[i:j]
                        if substring in substrings:
                            start_index = sequence.find(substring, substrings[substring] + 1)
                            if start_index != -1:
                                for k in range(start_index, start_index + len(substring)):
                                    result[k] = ""
                        else:
                            substrings[substring] = i
            
                return "".join(result)

            def get_vdj_region_column_name(df):
                """Returns the correct column name for VDJ region."""
                if 'V.D.J.REGION' in df.columns:
                    return 'V.D.J.REGION'
                elif 'V-D-J-REGION' in df.columns:
                    return 'V-D-J-REGION'
                else:
                    raise KeyError("Neither 'V.D.J.REGION' nor 'V-D-J-REGION' is present in the DataFrame.")


            # Function to process and save completed sequences
            def process_and_save_sequences(input_df, fasta_sequences, output_filename):
                """Processes input DataFrame and saves the completed sequences."""
            
                # Get the correct column name for VDJ region
                vdj_region_col = get_vdj_region_column_name(input_df)
                
                completed_sequences = []
                
                for _, row in input_df.iterrows():
                    query_sequence = row[vdj_region_col]
                    _, best_reference = find_best_alignment(fasta_sequences, query_sequence)
                
                    completed_sequence = extract_and_concatenate(query_sequence, best_reference) if best_reference else query_sequence
                    completed_sequences.append(completed_sequence)
                
                # Add completed sequences and remove repeats
                input_df["Completed_VH_Sequences"] = completed_sequences
                input_df["cleaned_sequence"] = input_df["Completed_VH_Sequences"].apply(keep_one_repeat)
                
                # Get the path to the Downloads folder and save the output CSV
                downloads_folder = str(Path.home() / "Downloads")
                output_csv = os.path.join(downloads_folder, output_filename)
                input_df.to_csv(output_csv, index=False)
                print(f"Output saved to {output_csv}")

            
            # Main function to process different groups
            def main():
                fasta_file_path = self.file_paths[8]  # Replace with actual path
                fasta_sequences = load_imgt_fasta(fasta_file_path)
            
                 # Process IgG1 clonally related to low mutated divergent shared IgE
                process_and_save_sequences(
                    clonally_related_low_mutated_shared_IgG1_seqid_in_IgG1_AA.copy(),
                    fasta_sequences,
                    "output_IgG1_clonally_related_to_low_mutated_divergent_shared_IgE_seqid_in_AA.csv"
                )
                
                # Process IgG1 clonally related to moderate mutated divergent shared IgE
                process_and_save_sequences(
                    clonally_related_moderate_mutated_shared_IgG1_seqid_in_IgG1_AA.copy(),
                    fasta_sequences,
                    "output_IgG1_clonally_related_to_moderate_mutated_divergent_shared_IgE_seqid_in_AA.csv"
                )
                    
                # Process IgG1 clonally related to high mutated divergent shared IgE
                process_and_save_sequences(
                    clonally_related_high_mutated_shared_IgG1_seqid_in_IgG1_AA.copy(),
                    fasta_sequences,
                    "output_IgG1_clonally_related_to_high_mutated_divergent_shared_IgE_seqid_in_AA.csv"
                )
            
                # Process top ten low mutated shared IgG1
                process_and_save_sequences(
                    stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten_tab.copy(),
                    fasta_sequences,
                    "output_stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten.csv"
                )
                
                 # Process top ten moderate mutated shared IgG1
                process_and_save_sequences(
                    stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten_tab.copy(),
                    fasta_sequences,
                    "output_stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten.csv"
                )
                
                # Process top ten high mutated shared IgG1
                process_and_save_sequences(
                    stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten_tab.copy(),
                    fasta_sequences,
                    "output_stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten.csv"
                )
                    
            if __name__ == "__main__":
                main()
            
                
            #Analyse the aminoacid composition from the shared IgE
            import re
            import pandas as pd
            import os
            from collections import Counter
            
            def extract_and_process_aa(sequence):
                if not isinstance(sequence, str):
                    return '', ''
                
                segments = [seg.strip() for seg in sequence.split('|') if seg.strip()]
                
                def process_mutation(mut_str):
                    """Process mutation string, treating repeated AAs as non-silent (A→A)"""
                    # Case 1: Explicit AA change (T105>S)
                    explicit_change = re.search(r'([A-Z])\d+>([A-Z])', mut_str)
                    if explicit_change:
                        return explicit_change.group(1), explicit_change.group(2)
                    
                    # Case 2: Repeated AA (A105; A105) - treated as non-silent (A→A)
                    repeated_aa = re.search(r'([A-Z])\d+;\s*([A-Z])\d+', mut_str)
                    if repeated_aa and repeated_aa.group(1) == repeated_aa.group(2):
                        return repeated_aa.group(1), repeated_aa.group(1)  # A→A
                    
                    return None, None
                
                germline = []
                mutated = []
                
                for segment in segments:
                    for mut_notation in [x.strip() for x in segment.split(',') if x.strip()]:
                        germ_aa, mut_aa = process_mutation(mut_notation)
                        if germ_aa is not None and mut_aa is not None:
                            germline.append(germ_aa)
                            mutated.append(mut_aa)
                
                return ''.join(germline), ''.join(mutated)
            
            # Create a set of unique seqid from the detected column
            mut_shared_IgE = stats_tab_mutated_shared_IgE_VDJ.copy()
            shared_IgE_VDJ_extract = set(mut_shared_IgE['seqid'])
            
            # Create a copy of the data frame
            IgE_AA_change_table1 = stats_IgE_Aminotab_change.copy()
            
            # Check for the correct column name in stats_IgE_Aminotab
            if 'Sequence.ID' in IgE_AA_change_table1.columns:
                column_seqid = 'Sequence.ID'
            elif 'Sequence ID' in IgE_AA_change_table1.columns:
                column_seqid = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present in the DataFrame.")
            
            # Update the comparison IgE DataFrame with boolean masks
            IgE_AA_change_table1['Match_in_mut_shared_IgE_seqid'] = IgE_AA_change_table1[column_seqid].isin(shared_IgE_VDJ_extract)
            
            # Ensure 'Match_in_mutated shared IgE' is a string and remove leading/trailing spaces
            IgE_AA_change_table1['Match_in_mut_shared_IgE_seqid'] = IgE_AA_change_table1['Match_in_mut_shared_IgE_seqid'].astype(str).str.strip()
            
            # Now filter the mutated shared IgE Aminotab change table
            IgE_AA_change_table2 = IgE_AA_change_table1.loc[IgE_AA_change_table1['Match_in_mut_shared_IgE_seqid'] == 'True'].copy()

            
            # Check for the correct column names in IgE_AA_change_table2
            if 'CDR3.IMGT' in IgE_AA_change_table2.columns:
                column_cdr3 = 'CDR3.IMGT'
            elif 'CDR3-IMGT' in IgE_AA_change_table2.columns:
                column_cdr3 = 'CDR3-IMGT'
            else:
                raise KeyError("Neither 'CDR3.IMGT' nor 'CDR3-IMGT' column is present in the DataFrame.")
            
            if 'CDR2.IMGT' in IgE_AA_change_table2.columns:
                column_cdr2 = 'CDR2.IMGT'
            elif 'CDR2-IMGT' in IgE_AA_change_table2.columns:
                column_cdr2 = 'CDR2-IMGT'
            else:
                raise KeyError("Neither 'CDR2.IMGT' nor 'CDR2-IMGT' column is present in the DataFrame.")
            
            # Lists to store the germline and mutated sequences
            germline_column_cdr3 = []
            mutated_column_cdr3 = []
            germline_column_cdr2 = []
            mutated_column_cdr2 = []
            
            # Apply the function to both CDR3 and CDR2 columns
            for cdr3, cdr2 in zip(IgE_AA_change_table2[column_cdr3], IgE_AA_change_table2[column_cdr2]):
                # Process CDR3
                germline_cdr3, mutated_cdr3 = extract_and_process_aa(cdr3)
                germline_column_cdr3.append(germline_cdr3)
                mutated_column_cdr3.append(mutated_cdr3)
            
                # Process CDR2
                germline_cdr2, mutated_cdr2 = extract_and_process_aa(cdr2)
                germline_column_cdr2.append(germline_cdr2)
                mutated_column_cdr2.append(mutated_cdr2)
            
            # Add the extracted data to the DataFrame
            IgE_AA_change_table2.loc[:, 'Germline_CDR3'] = germline_column_cdr3
            IgE_AA_change_table2.loc[:, 'Mutated_CDR3'] = mutated_column_cdr3
            IgE_AA_change_table2.loc[:, 'Germline_CDR2'] = germline_column_cdr2
            IgE_AA_change_table2.loc[:, 'Mutated_CDR2'] = mutated_column_cdr2
                        
            # Concatenate the amino acid sequences for each region into a single string (entire column)
            Germline_CDR3_concat = ''.join(IgE_AA_change_table2['Germline_CDR3'])
            Mutated_CDR3_concat = ''.join(IgE_AA_change_table2['Mutated_CDR3'])
            Germline_CDR2_concat = ''.join(IgE_AA_change_table2['Germline_CDR2'])
            Mutated_CDR2_concat = ''.join(IgE_AA_change_table2['Mutated_CDR2'])
            
            # Define EMBOSS amino acid categories
            emboss_categories = {
                "Tiny": set("ACGST"),
                "Small": set("ABCDGNPSTV"),
                "Aliphatic": set("AILV"),
                "Aromatic": set("FHWY"),
                "Non-polar": set("ACFGILMPVWY"),
                "Polar": set("DEHKNQRSTZ"),
                "Charged": set("BDEHKRZ"),
                "Basic": set("KRH"),
                "Acidic": set("BDEZ")
            }
            
            # Calculate amino acid composition by EMBOSS categories
            def calculate_emboss_composition(sequence):
                if not sequence:
                    return {category: 0 for category in emboss_categories.keys()}
                
                # Count ALL characters for total length (including non-standard AAs)
                total_aa = len(sequence)
                                          
                # Define standard amino acids (20 standard + B/Z for ambiguous)
                standard_aas = set("ACDEFGHIKLMNPQRSTVWYBZ")
            
                # Initialize counters
                aa_count = Counter()
                category_counts = {category: 0 for category in emboss_categories.keys()}
                
            
                # Process each amino acid
                for aa in sequence:
                    # Count all AAs (including non-standard)
                    aa_count[aa] += 1
                    
                    # Only categorize standard AAs
                    if aa in standard_aas:
                        for category, aa_set in emboss_categories.items():
                            if aa in aa_set:
                                category_counts[category] += 1
                
                # Calculate percentages (relative to TOTAL sequence length)
                category_percentages = {
                    category: (count / total_aa) * 100 
                    for category, count in category_counts.items()
                }
                
                return category_percentages
            
            # Calculate the composition for each region
            Germline_CDR3_composition = calculate_emboss_composition(Germline_CDR3_concat)
            Mutated_CDR3_composition = calculate_emboss_composition(Mutated_CDR3_concat)
            Germline_CDR2_composition = calculate_emboss_composition(Germline_CDR2_concat)
            Mutated_CDR2_composition = calculate_emboss_composition(Mutated_CDR2_concat)
            
            # Combine all the compositions into a dictionary
            composition_results = {
                'Germline_CDR3': Germline_CDR3_composition,
                'Mutated_CDR3': Mutated_CDR3_composition,
                'Germline_CDR2': Germline_CDR2_composition,
                'Mutated_CDR2': Mutated_CDR2_composition,
            }
            
            # Dynamically determine the Downloads folder
            def get_downloads_folder():
                if os.name == "nt":  # Windows
                    return os.path.join(os.environ["USERPROFILE"], "Downloads")
                else:  # macOS/Linux
                    return os.path.join(os.path.expanduser("~"), "Downloads")
            
            downloads_folder = get_downloads_folder()
            output_file_path = os.path.join(downloads_folder, "Shared_IgE_AminoAcid_Composition.xlsx")
            
            # Convert composition results to DataFrame
            composition_df = pd.DataFrame(composition_results)
            
            # Save the DataFrame as an Excel file
            composition_df.to_excel(output_file_path, index=True)
            print(f"File saved successfully in {output_file_path}")

    
            #Analyse the aminoacid composition from the shared IgG1
            import re
            import pandas as pd
            import os
            from collections import Counter
            
            # Define the function to extract and process the amino acid sequences
            def extract_and_process_aa(sequence):
                if not isinstance(sequence, str):
                    return '', ''
                
                segments = [seg.strip() for seg in sequence.split('|') if seg.strip()]
                
                def process_mutation(mut_str):
                    """Process mutation string, treating repeated AAs as non-silent (A→A)"""
                    # Case 1: Explicit AA change (T105>S)
                    explicit_change = re.search(r'([A-Z])\d+>([A-Z])', mut_str)
                    if explicit_change:
                        return explicit_change.group(1), explicit_change.group(2)
                    
                    # Case 2: Repeated AA (A105; A105) - treated as non-silent (A→A)
                    repeated_aa = re.search(r'([A-Z])\d+;\s*([A-Z])\d+', mut_str)
                    if repeated_aa and repeated_aa.group(1) == repeated_aa.group(2):
                        return repeated_aa.group(1), repeated_aa.group(1)  # A→A
                    
                    return None, None
                
                germline = []
                mutated = []
                
                for segment in segments:
                    for mut_notation in [x.strip() for x in segment.split(',') if x.strip()]:
                        germ_aa, mut_aa = process_mutation(mut_notation)
                        if germ_aa is not None and mut_aa is not None:
                            germline.append(germ_aa)
                            mutated.append(mut_aa)
                
                return ''.join(germline), ''.join(mutated)
            
            # Create a set of unique seqid from the detected column
            mut_shared_IgG1 = stats_tab_mutated_shared_IgG1_VDJ.copy()
            shared_IgG1_VDJ_extract = set(mut_shared_IgG1['seqid'])
            
            # Create a copy of the data frame
            IgG1_AA_change_table1 = stats_IgG1_Aminotab_change.copy()
            
            # Check for the correct column name in stats_IgG1_Aminotab
            if 'Sequence.ID' in IgG1_AA_change_table1.columns:
                column_seqid = 'Sequence.ID'
            elif 'Sequence ID' in IgG1_AA_change_table1.columns:
                column_seqid = 'Sequence ID'
            else:
                raise KeyError("Neither 'Sequence.ID' nor 'Sequence ID' column is present in the DataFrame.")
            
            # Update the comparison IgG1 DataFrame with boolean masks
            IgG1_AA_change_table1['Match_in_mut_shared_IgG1_seqid'] = IgG1_AA_change_table1[column_seqid].isin(shared_IgG1_VDJ_extract)
            
            # Ensure 'Match_in_mutated shared IgG1' is a string and remove leading/trailing spaces
            IgG1_AA_change_table1['Match_in_mut_shared_IgG1_seqid'] = IgG1_AA_change_table1['Match_in_mut_shared_IgG1_seqid'].astype(str).str.strip()
            
            # Now filter the mutated shared IgG1 Aminotab change table
            IgG1_AA_change_table2 = IgG1_AA_change_table1.loc[IgG1_AA_change_table1['Match_in_mut_shared_IgG1_seqid'] == 'True'].copy()
            
            # Check for the correct column names in IgG1_AA_change_table2
            if 'CDR3.IMGT' in IgG1_AA_change_table2.columns:
                column_cdr3 = 'CDR3.IMGT'
            elif 'CDR3-IMGT' in IgG1_AA_change_table2.columns:
                column_cdr3 = 'CDR3-IMGT'
            else:
                raise KeyError("Neither 'CDR3.IMGT' nor 'CDR3-IMGT' column is present in the DataFrame.")
            
            if 'CDR2.IMGT' in IgG1_AA_change_table2.columns:
                column_cdr2 = 'CDR2.IMGT'
            elif 'CDR2-IMGT' in IgG1_AA_change_table2.columns:
                column_cdr2 = 'CDR2-IMGT'
            else:
                raise KeyError("Neither 'CDR2.IMGT' nor 'CDR2-IMGT' column is present in the DataFrame.")
            
            # Lists to store the germline and mutated sequences
            germline_column_cdr3 = []
            mutated_column_cdr3 = []
            germline_column_cdr2 = []
            mutated_column_cdr2 = []
            
            # Apply the function to both CDR3 and CDR2 columns
            for cdr3, cdr2 in zip(IgG1_AA_change_table2[column_cdr3], IgG1_AA_change_table2[column_cdr2]):
                # Process CDR3
                germline_cdr3, mutated_cdr3 = extract_and_process_aa(cdr3)
                germline_column_cdr3.append(germline_cdr3)
                mutated_column_cdr3.append(mutated_cdr3)
            
                # Process CDR2
                germline_cdr2, mutated_cdr2 = extract_and_process_aa(cdr2)
                germline_column_cdr2.append(germline_cdr2)
                mutated_column_cdr2.append(mutated_cdr2)
            
            # Add the extracted data to the DataFrame
            IgG1_AA_change_table2.loc[:, 'Germline_CDR3'] = germline_column_cdr3
            IgG1_AA_change_table2.loc[:, 'Mutated_CDR3'] = mutated_column_cdr3
            IgG1_AA_change_table2.loc[:, 'Germline_CDR2'] = germline_column_cdr2
            IgG1_AA_change_table2.loc[:, 'Mutated_CDR2'] = mutated_column_cdr2

            # Concatenate the amino acid sequences for each region into a single string (entire column)
            Germline_CDR3_concat = ''.join(IgG1_AA_change_table2['Germline_CDR3'])
            Mutated_CDR3_concat = ''.join(IgG1_AA_change_table2['Mutated_CDR3'])
            Germline_CDR2_concat = ''.join(IgG1_AA_change_table2['Germline_CDR2'])
            Mutated_CDR2_concat = ''.join(IgG1_AA_change_table2['Mutated_CDR2'])
            
            # Define EMBOSS amino acid categories
            emboss_categories = {
                "Tiny": set("ACGST"),
                "Small": set("ABCDGNPSTV"),
                "Aliphatic": set("AILV"),
                "Aromatic": set("FHWY"),
                "Non-polar": set("ACFGILMPVWY"),
                "Polar": set("DEHKNQRSTZ"),
                "Charged": set("BDEHKRZ"),
                "Basic": set("KRH"),
                "Acidic": set("BDEZ")
            }
            
            # Calculate amino acid composition by EMBOSS categories
            def calculate_emboss_composition(sequence):
                if not sequence:
                    return {category: 0 for category in emboss_categories.keys()}
                
                # Count ALL characters for total length (including non-standard AAs)
                total_aa = len(sequence)
                                      
                # Define standard amino acids (20 standard + B/Z for ambiguous)
                standard_aas = set("ACDEFGHIKLMNPQRSTVWYBZ")
            
                # Initialize counters
                aa_count = Counter()
                category_counts = {category: 0 for category in emboss_categories.keys()}
                
            
                # Process each amino acid
                for aa in sequence:
                    # Count all AAs (including non-standard)
                    aa_count[aa] += 1
                    
                    # Only categorize standard AAs
                    if aa in standard_aas:
                        for category, aa_set in emboss_categories.items():
                            if aa in aa_set:
                                category_counts[category] += 1
                
                # Calculate percentages (relative to TOTAL sequence length)
                category_percentages = {
                    category: (count / total_aa) * 100 
                    for category, count in category_counts.items()
                }
                
                return category_percentages
            
            # Calculate the composition for each region
            Germline_CDR3_composition = calculate_emboss_composition(Germline_CDR3_concat)
            Mutated_CDR3_composition = calculate_emboss_composition(Mutated_CDR3_concat)
            Germline_CDR2_composition = calculate_emboss_composition(Germline_CDR2_concat)
            Mutated_CDR2_composition = calculate_emboss_composition(Mutated_CDR2_concat)
            
            # Combine all the compositions into a dictionary
            composition_results = {
                'Germline_CDR3': Germline_CDR3_composition,
                'Mutated_CDR3': Mutated_CDR3_composition,
                'Germline_CDR2': Germline_CDR2_composition,
                'Mutated_CDR2': Mutated_CDR2_composition,
            }
            
            # Dynamically determine the Downloads folder
            def get_downloads_folder():
                if os.name == "nt":  # Windows
                    return os.path.join(os.environ["USERPROFILE"], "Downloads")
                else:  # macOS/Linux
                    return os.path.join(os.path.expanduser("~"), "Downloads")
            
            downloads_folder = get_downloads_folder()
            output_file_path = os.path.join(downloads_folder, "Shared_IgG1_AminoAcid_Composition.xlsx")
            
            # Convert composition results to DataFrame
            composition_df = pd.DataFrame(composition_results)
            
            # Save the DataFrame as an Excel file
            composition_df.to_excel(output_file_path, index=True)
            print(f"File saved successfully in {output_file_path}")
                    
            
            # Create a custom function to capture print statements
            class PrintCapture:
                def __init__(self):
                    self.output = StringIO()  # Create an in-memory text stream
                    self.saved_stdout = sys.stdout  # Save the original stdout
                    sys.stdout = self.output  # Redirect stdout to the in-memory text stream
            
                def get_output(self):
                    return self.output.getvalue()  # Return the captured output
            
                def reset(self):
                    sys.stdout = self.saved_stdout  # Restore the original stdout
            
            # Function to perform the operations and capture print output
            def your_script_function():
                # Example print statements (replace these with your actual outputs)
                print(f'Number of IgE clones: {unique_strings_IgE_productive_VDJ}')
                print(f'Number of IgG1 clones: {unique_strings_IgG1_productive_VDJ}')
                print(f'Number of shared IgE clones: {unique_strings_shared_IgE}')
                print(f'percentage of shared IgE clones: {percentage_shared_IgE_clones}')
                print(f'Number of unique IgE clones: {unique_strings_unique_IgE}')
                print(f'percentage of unique IgE clones: {percentage_unique_IgE_clones}')
                print(f'Number of shared IgG1 clones: {unique_strings_shared_IgG1}')
                print(f'percentage of shared IgG1 clones: {percentage_shared_IgG1_clones}')
                print(f'Number of unique IgG1 clones: {unique_strings_unique_IgG1}')
                print(f'percentage of unique IgG1 clones: {percentage_unique_IgG1_clones}')
                print(f'Mutated IgE_biased_average percentage of CDR2 non silent mutation: {mean_CDR2_stats_tab_mutated_IgE_bias}')
                print(f'Mutated IgE_biased_average percentage of CDR3 non silent mutation: {mean_CDR3_stats_tab_mutated_IgE_bias}')
                print(f'Mutated IgG1_biased_average percentage of CDR2 non silent mutation: {mean_CDR2_stats_tab_mutated_IgG1_bias}')
                print(f'Mutated IgG1_biased_average percentage of CDR3 non silent mutation: {mean_CDR3_stats_tab_mutated_IgG1_bias}')
                print(f'Mutated unique_IgE_average percentage of CDR2 non silent mutation: {mean_CDR2_stats_tab_mutated_Unique_IgE}')
                print(f'Mutated unique_IgE__average percentage of CDR3 non silent mutation: {mean_CDR3_stats_tab_mutated_Unique_IgE}')
                print(f'Unmutated IgE_biased_average percentage of CDR2 non silent mutation: {mean_CDR2_stats_tab_unmutated_IgE_bias}')
                print(f'Unmutated IgE_biased_average percentage of CDR3 non silent mutation: {mean_CDR3_stats_tab_unmutated_IgE_bias}')
                print(f'Unmutated IgG1_biased_average percentage of CDR2 non silent mutation: {mean_CDR2_stats_tab_unmutated_IgG1_bias}')
                print(f'Unmutated IgG1_biased_average percentage of CDR3 non silent mutation: {mean_CDR3_stats_tab_unmutated_IgG1_bias}')
                print(f'Unmutated unique_IgE_average percentage of CDR2 non silent mutation: {mean_CDR2_stats_tab_unmutated_Unique_IgE}')
                print(f'Unmutated unique_IgE_average percentage of CDR3 non silent mutation: {mean_CDR3_stats_tab_unmutated_Unique_IgE}')
                print(f'Number of mutated IgE-biased clones: {unique_strings_mutated_IgE_bias_VDJ}')
                print(f'percentage of mutated IgE-biased clones: {percentage_mutated_IgE_bias_clones}')
                print(f'mean copy number of mutated IgE-biased clones: {mean_copynumber_mutated_IgE_bias_VDJ}')
                print(f'Number of mutated IgG1-biased clones: {unique_strings_mutated_IgG1_bias_VDJ}')
                print(f'percentage of mutated IgG1-biased clones: {percentage_mutated_IgG1_bias_clones}')
                print(f'mean copy number of mutated IgG1-biased clones: {mean_copynumber_mutated_IgG1_bias_VDJ}')
                print(f'Number of mutated unique IgE clones: {unique_strings_mutated_unique_IgE_VDJ}')
                print(f'percentage of mutated unique IgE clones: {percentage_mutated_unique_IgE_clones}')
                print(f'mean copy number of mutated unique IgE clones: {mean_copynumber_mutated_unique_IgE_VDJ}')
                print(f'Number of mutated Total IgE clones: {unique_strings_mutated_Total_IgE_VDJ}')
                print(f'percentage of mutated Total IgE clones: {percentage_mutated_Total_IgE_clones}')
                print(f'mean copy number of mutated Total IgE clones: {mean_copynumber_mutated_Total_IgE_VDJ}')
                print(f'Number of unmutated IgE-biased clones: {unique_strings_unmutated_IgE_bias_VDJ}')
                print(f'percentage of unmutated IgE-biased clones: {percentage_unmutated_IgE_bias_clones}')
                print(f'mean copy number of unmutated IgE-biased clones: {mean_copynumber_unmutated_IgE_bias_VDJ}')
                print(f'Number of unmutated IgG1-biased clones: {unique_strings_unmutated_IgG1_bias_VDJ}')
                print(f'percentage of unmutated IgG1-biased clones: {percentage_unmutated_IgG1_bias_clones}')
                print(f'mean copy number of unmutated IgG1-biased clones: {mean_copynumber_unmutated_IgG1_bias_VDJ}')
                print(f'Number of unmutated unique IgE clones: {unique_strings_unmutated_unique_IgE_VDJ}')
                print(f'percentage of unmutated unique IgE clones: {percentage_unmutated_unique_IgE_clones}')
                print(f'mean copy number of unmutated unique IgE clones: {mean_copynumber_unmutated_unique_IgE_VDJ}')
                print(f'Number of unmutated Total IgE clones: {unique_strings_unmutated_Total_IgE_VDJ}')
                print(f'percentage of unmutated Total IgE clones: {percentage_unmutated_Total_IgE_clones}')
                print(f'mean copy number of unmutated Total IgE clones: {mean_copynumber_unmutated_Total_IgE_VDJ}')
                print(f'Number of mutated Total IgG1 clones: {unique_strings_mutated_Total_IgG1_VDJ}')
                print(f'Number of unmutated Total IgG1 clones: {unique_strings_unmutated_Total_IgG1_VDJ}')
                print(f'Percentage of mutated Total IgG1 clones: {percentage_mutated_Total_IgG1_clones:.2f}%')
                print(f'Mean copy number of mutated Total IgG1 clones: {mean_copynumber_mutated_Total_IgG1_VDJ}')
                print(f'Percentage of unmutated Total IgG1 clones: {percentage_unmutated_Total_IgG1_clones:.2f}%')
                print(f'Mean copy number of unmutated Total IgG1 clones: {mean_copynumber_unmutated_Total_IgG1_VDJ}')
                print(f'mean copynumber low mutated IgE: {mean1}')
                print(f'mean copynumber low mutated IgG1: {mean2}')
                print(f'mean copynumber moderate mutated IgE: {mean3}')
                print(f'mean copynumber moderate mutated IgG1: {mean4}')
                print(f'mean copynumber high mutated IgE: {mean5}')
                print(f'mean copynumber high mutated IgG1: {mean6}')
                print(f'Mean CDR3 netcharge of low mutated IgE: {mean_stats_tab_low_mutated_shared_IgE_VDJ_sorted_topten}')
                print(f'Mean CDR3 netcharge of moderate mutated IgE: {mean_stats_tab_moderate_mutated_shared_IgE_VDJ_sorted_topten}')
                print(f'Mean CDR3 netcharge of high mutated IgE: {mean_stats_tab_high_mutated_shared_IgE_VDJ_sorted_topten}')
                print(f'Mean CDR3 netcharge of low mutated IgG1: {mean_stats_tab_low_mutated_shared_IgG1_VDJ_sorted_topten}')
                print(f'Mean CDR3 netcharge of moderate mutated IgG1: {mean_stats_tab_moderate_mutated_shared_IgG1_VDJ_sorted_topten}')
                print(f'Mean CDR3 netcharge of high mutated IgG1: {mean_stats_tab_high_mutated_shared_IgG1_VDJ_sorted_topten}')
                print(f'low mutated IgG1/IgE netcharge ratio: {ratio_stats_tab_low_mutated_shared_IgG1_IgE_VDJ_sorted_topten}')
                print(f'moderate mutated IgG1/IgE netcharge ratio: {ratio_stats_tab_moderate_mutated_shared_IgG1_IgE_VDJ_sorted_topten}')
                print(f'high mutated IgG1/IgE netcharge ratio: {ratio_stats_tab_high_mutated_shared_IgG1_IgE_VDJ_sorted_topten}')
                print(f'Number of divergent CDR3 low mutated IgE clones: {unique_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone}')
                print(f'Number of divergent CDR3 moderate mutated IgE clones: {unique_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone}')
                print(f'Number of divergent CDR3 high mutated IgE clones: {unique_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone}')
                print(f'Number of non_divergent CDR3 low mutated IgE clones: {unique_non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone}')
                print(f'Number of non_divergent CDR3 moderate mutated IgE clones: {unique_non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone}')
                print(f'Number of non_divergent CDR3 high mutated IgE clones: {unique_non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone}')
                print(f'mean copynumber of divergent CDR3 low mutated IgE: {mean_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone}')
                print(f'mean copynumber of divergent CDR3 moderate mutated IgE: {mean_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone}')
                print(f'mean copynumber of divergent CDR3 high mutated IgE: {mean_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone}')
                print(f'mean copynumber of non_divergent CDR3 low mutated IgE: {mean_non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone}')
                print(f'mean copynumber of non_divergent CDR3 moderate mutated IgE: {mean_non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone}')
                print(f'mean copynumber of non_divergent CDR3 high mutated IgE: {mean_non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone}')
                print(f'sum copynumber of divergent CDR3 low mutated IgE: {sum_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone}')
                print(f'sum copynumber of divergent CDR3 moderate mutated IgE: {sum_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone}')
                print(f'sum copynumber of divergent CDR3 high mutated IgE: {sum_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone}')
                print(f'sum copynumber of non_divergent CDR3 low mutated IgE: {sum_non_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3mut_clone}')
                print(f'sum copynumber of non_divergent CDR3 moderate mutated IgE: {sum_non_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3mut_clone}')
                print(f'sum copynumber of non_divergent CDR3 high mutated IgE: {sum_non_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3mut_clone}')
                print(f'mean copynumber of IgG1 clonally related to low mutated IgE: {mean_clonally_related_low_mutated_shared_IgG1}')
                print(f'mean copynumber of IgG1 clonally related to moderate mutated IgE: {mean_clonally_related_moderate_mutated_shared_IgG1}')
                print(f'mean copynumber of IgG1 clonally related to high mutated IgE: {mean_clonally_related_high_mutated_shared_IgG1}')
                print(f'sum copynumber of IgG1 clonally related to low mutated IgE: {sum_clonally_related_low_mutated_shared_IgG1}')
                print(f'sum copynumber of IgG1 clonally related to moderate mutated IgE: {sum_clonally_related_moderate_mutated_shared_IgG1}')
                print(f'sum copynumber of IgG1 clonally related to high mutated IgE: {sum_clonally_related_high_mutated_shared_IgG1}')
                print(f'Number of divergent CDR3CDR2 low mutated IgE clones: {unique_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone}')
                print(f'Number of divergent CDR3CDR2 moderate mutated IgE clones: {unique_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone}')
                print(f'Number of divergent CDR3CDR2 high mutated IgE clones: {unique_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone}')
                print(f'mean copynumber of divergent CDR3CDR2 low mutated IgE: {mean_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone}')
                print(f'mean copynumber of divergent CDR3CDR2 moderate mutated IgE: {mean_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone}')
                print(f'mean copynumber of divergent CDR3CDR2 high mutated IgE: {mean_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone}')
                print(f'sum copynumber of divergent CDR3CDR2 low mutated IgE: {sum_divergent_tab_low_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone}')
                print(f'sum copynumber of divergent CDR3CDR2 moderate mutated IgE: {sum_divergent_tab_moderate_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone}')
                print(f'sum copynumber of divergent CDR3CDR2 high mutated IgE: {sum_divergent_tab_high_mutated_shared_IgE_VDJ_CDR3CDR2mut_clone}')
                print(f'mean copynumber of IgG1 clonally  related to low mutated IgE (CDR3CDR2): {mean_clonally_related_low_mutated_shared_CDR3CDR2_IgG1}')
                print(f'mean copynumber of IgG1 clonally related to moderate mutated IgE (CDR3CDR2): {mean_clonally_related_moderate_mutated_shared_CDR3CDR2_IgG1}')
                print(f'mean copynumber of IgG1 clonally related to high mutated IgE (CDR3CDR2): {mean_clonally_related_high_mutated_shared_CDR3CDR2_IgG1}')
            
            # Create an instance of PrintCapture to capture print output
            capture = PrintCapture()
            
            # Run the function that generates print outputs
            your_script_function()
            
            # Get all the captured print outputs
            captured_output = capture.get_output()
            
            # Reset the capture (restore stdout)
            capture.reset()
            
            # Now you can save the captured output to an Excel file
            output_lines = captured_output.splitlines()  # Split the captured output into lines
            
            # Create a DataFrame from the captured output
            df = pd.DataFrame(output_lines, columns=["Captured Output"])
            
            # Determine the Downloads folder of the current user
            downloads_folder = os.path.join(os.path.expanduser("~"), "Downloads")
            output_excel_file = os.path.join(downloads_folder, "captured_output.xlsx")
            
            # Create a DataFrame from the captured output
            df = pd.DataFrame(output_lines, columns=["Captured Output"])
            
            # Determine the Downloads folder of the current user
            downloads_folder = os.path.join(os.path.expanduser("~"), "Downloads")
            output_excel_file = os.path.join(downloads_folder, "captured_output.xlsx")
            
            # Save the DataFrame to an Excel file in the Downloads folder
            df.to_excel(output_excel_file, index=False)
            
            print(f"Captured output has been saved to {output_excel_file}")
            
            self.analysis_completed.emit("Analysis completed successfully!")
        except Exception as e:
            self.error_occurred.emit(str(e))


class Screen2(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Plots Display")
        self.setGeometry(200, 200, 800, 600)
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.current_index = 0
        self.figures = []  # To store all active figures
        self.canvas = None

        self.load_all_plots()

        # Next Plot Button
        next_button = QPushButton("Next Plot")
        next_button.clicked.connect(self.show_next_plot)
        self.layout.addWidget(next_button)

    def load_all_plots(self):
        """Load all active matplotlib figures."""
        for fignum in plt.get_fignums():  # Get all active figure numbers
            fig = plt.figure(fignum)      # Get the figure object by number
            self.figures.append(fig)

        if self.figures:
            self.display_plot(0)  # Display the first plot initially

    def display_plot(self, index):
        """Display the plot at the given index."""
        if self.canvas:  # Remove existing canvas if present
            self.layout.removeWidget(self.canvas)
            self.canvas.deleteLater()

        self.canvas = FigureCanvas(self.figures[index])
        self.layout.insertWidget(0, self.canvas)

    def show_next_plot(self):
        """Show the next plot in the list."""
        if not self.figures:
            return
        self.current_index = (self.current_index + 1) % len(self.figures)
        self.display_plot(self.current_index)



class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle(" BCR Analysis Tool")
        self.setGeometry(100, 100, 1000, 600)
        
        # Create main layout
        main_layout = QVBoxLayout()
        form_layout = QVBoxLayout()

        # Apply black, grey, and white palette with gradient
        palette = QPalette()
        gradient = QLinearGradient(0, 0, 1, 1)
        gradient.setColorAt(0.0, QColor(10, 10, 10))
        gradient.setColorAt(1.0, QColor(40, 40, 40))
        palette.setBrush(QPalette.Window, QBrush(gradient))
        palette.setColor(QPalette.WindowText, QColor(255, 255, 255))
        palette.setColor(QPalette.Base, QColor(20, 20, 20))
        palette.setColor(QPalette.Text, QColor(200, 200, 200))
        self.setPalette(palette)

        # Header with glow effect
        header_label = QLabel("CompIgS")
        header_label.setAlignment(Qt.AlignCenter)
        header_label.setFont(QFont("Orbitron", 22, QFont.Bold))
        header_label.setStyleSheet("""
            QLabel {
                color: #cccccc;
                text-shadow: 2px 2px 4px rgba(100, 100, 100, 0.5);
            }
        """)
        main_layout.addWidget(header_label)

        # File input fields
        self.file_inputs = []
        file_labels = [
            "(Ig1)IgE IMGT StatClonotype Output File:",
            "(Ig2)IgG1 IMGT StatClonotype Output File:",
            "(Ig1)IgE 8-V-region-nt-mutation-statistics table:",
            "(Ig2)IgG1 8-V-region-nt-mutation-statistics table:",
            "(Ig1)IgE 5-AA-sequences table:",
            "(Ig2)IgG1 5-AA-sequences table:",
            "(Ig1)IgE 7-V-REGION-mutation-and-AA-change-table:",
            "(Ig2)IgG1 7-V-REGION-mutation-and-AA-change-table:",
            "IMGT mouse VH reference:"
        ]

        for label_text in file_labels:
            row_layout = QHBoxLayout()
            label = QLabel(label_text)
            label.setFont(QFont("Courier", 12))
            label.setStyleSheet("color: #999999;")
            line_edit = QLineEdit()
            line_edit.setStyleSheet("""
                QLineEdit {
                    background-color: rgba(20, 20, 20, 0.9);
                    color: #cccccc;
                    border: 1px solid #666666;
                    padding: 6px;
                    border-radius: 5px;
                }
                QLineEdit:hover {
                    border: 1px solid #999999;
                }
            """)
            browse_button = QPushButton("Browse")
            browse_button.setStyleSheet("""
                QPushButton {
                    background-color: rgba(100, 100, 100, 0.2);
                    color: #cccccc;
                    border: 1px solid #666666;
                    padding: 6px 12px;
                    border-radius: 5px;
                }
                QPushButton:hover {
                    background-color: rgba(150, 150, 150, 0.4);
                }
            """)
            browse_button.clicked.connect(lambda _, le=line_edit: self.browse_file(le))

            row_layout.addWidget(label)
            row_layout.addWidget(line_edit)
            row_layout.addWidget(browse_button)

            form_layout.addLayout(row_layout)
            self.file_inputs.append(line_edit)

        main_layout.addLayout(form_layout)

        # Horizontal separator
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        separator.setStyleSheet("background-color: #666666; height: 2px;")
        main_layout.addWidget(separator)

        # Progress bar with glow effect
        self.progress_bar = QProgressBar(self)
        self.progress_bar.setStyleSheet("""
            QProgressBar {
                border: 2px solid #666666;
                border-radius: 6px;
                text-align: center;
                background-color: rgba(20, 20, 20, 0.9);
                color: #cccccc;
            }
            QProgressBar::chunk {
                background-color: #666666;
                width: 20px;
            }
        """)

        self.status_label = QLabel("Status: Waiting for input")
        self.status_label.setAlignment(Qt.AlignCenter)
        self.status_label.setFont(QFont("Courier", 12))
        self.status_label.setStyleSheet("color: #999999;")

        main_layout.addWidget(self.progress_bar)
        main_layout.addWidget(self.status_label)

        # Buttons with monochrome style
        button_layout = QHBoxLayout()
        start_button = QPushButton("Start Analysis")
        start_button.setStyleSheet("""
            QPushButton {
                background-color: #666666;
                color: white;
                padding: 10px 20px;
                font-size: 14px;
                border-radius: 8px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #999999;
            }
        """)
        start_button.clicked.connect(self.start_analysis)

        self.goto_screen2_button = QPushButton("view Plots", self)
        self.goto_screen2_button.setStyleSheet("""
            QPushButton {
                background-color: #444444;
                color: white;
                padding: 10px 20px;
                font-size: 14px;
                border-radius: 8px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #777777;
            }
        """)
        self.goto_screen2_button.clicked.connect(self.gotoScreen2)

        # Clear All button
        clear_all_button = QPushButton("Clear All")
        clear_all_button.setStyleSheet("""
            QPushButton {
                background-color: #555555;
                color: white;
                padding: 10px 20px;
                font-size: 14px;
                border-radius: 8px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #888888;
            }
        """)
        clear_all_button.clicked.connect(self.clear_all_inputs)
        
        # Exit button
        exit_button = QPushButton("Exit")
        exit_button.setStyleSheet("""
            QPushButton {
                background-color: #333333;
                color: white;
                padding: 10px 20px;
                font-size: 14px;
                border-radius: 8px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #666666;
            }
        """)
        exit_button.clicked.connect(self.close)

        button_layout.addWidget(start_button)
        button_layout.addWidget(self.goto_screen2_button)
        button_layout.addWidget(exit_button)
        button_layout.addWidget(clear_all_button)
        main_layout.addLayout(button_layout)

        # Set the central widget
        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)
    def browse_file(self, line_edit):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File")
        if file_path:
            line_edit.setText(file_path)

    def start_analysis(self):
        file_paths = [le.text() for le in self.file_inputs]

        if all(file_paths):
            self.thread = AnalysisThread(file_paths)
            self.thread.progress_updated.connect(self.update_progress)
            self.thread.analysis_completed.connect(self.on_analysis_complete)
            self.thread.error_occurred.connect(self.on_error)

            self.thread.start()
            self.status_label.setText("Status: Analysis in progress...")
        else:
            QMessageBox.warning(self, "Input Error", "Please provide all required files.")
            self.status_label.setText("Error: Missing required files.")

    def update_progress(self, value):
        self.progress_bar.setValue(value)

    def on_analysis_complete(self, message):
        self.status_label.setText(message)
        QMessageBox.information(self, "Success", message)

    def on_error(self, error_message):
        self.status_label.setText(f"Error: {error_message}")
        QMessageBox.critical(self, "Error", error_message)
        
    def gotoScreen2(self):
        self.screen2 = Screen2()  # Error: Screen2 is not defined yet
        self.screen2.exec_()
        
    def clear_all_inputs(self):
        """Clears all input file paths from the input fields."""
        for line_edit in self.file_inputs:
            line_edit.clear()  # Clear the text in each QLineEdit
        plt.close('all')
        self.status_label.setText("Status: All input fields cleared.")
        
# Run the application
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())


# In[ ]:





# In[ ]:




