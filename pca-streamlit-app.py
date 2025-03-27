import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
import io
import base64
from matplotlib.backends.backend_pdf import PdfPages
import re

st.set_page_config(
    page_title="Gene Expression PCA Analysis Tool",
    page_icon="🧬",
    layout="wide"
)

# Custom CSS for styling
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #4682B4;
        text-align: center;
        margin-bottom: 1rem;
    }
    .sub-header {
        font-size: 1.5rem;
        color: #4682B4;
        margin-bottom: 1rem;
    }
    .info-box {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        margin-bottom: 1rem;
    }
    .interpretation {
        background-color: #e6f3ff;
        padding: 1rem;
        border-radius: 0.5rem;
        margin-top: 1rem;
    }
</style>
""", unsafe_allow_html=True)

# Header
st.markdown("<h1 class='main-header'>Gene Expression PCA Analysis Tool</h1>", unsafe_allow_html=True)

# Introduction
with st.expander("About this tool", expanded=False):
    st.markdown("""
    ### Purpose
    This tool performs Principal Component Analysis (PCA) on gene expression data to help identify patterns and relationships between samples.
    
    ### Features
    - Automatic sample group detection
    - Flexible data format support
    - High-resolution plot exports (PNG, JPG, PDF)
    - Comprehensive PCA interpretation
    - Sample and gene contribution analysis
    
    ### How to use
    1. Upload your gene expression data (CSV, TSV, or Excel)
    2. Configure sample grouping and visualization settings
    3. Analyze the PCA plot and interpretation
    4. Download high-resolution images for publication
    """)

# File upload
st.markdown("<h2 class='sub-header'>1. Upload your data</h2>", unsafe_allow_html=True)
uploaded_file = st.file_uploader("Upload gene expression data file", type=["csv", "tsv", "txt", "xlsx", "xls"])

# Demo data option
use_demo = st.checkbox("Use demo data (gene expression from sample image)")

# Process uploaded data
if uploaded_file is not None or use_demo:
    try:
        if use_demo:
            # Create demo data based on the image
            data = {
                'Gene Symbol': ['unc-54', 'unc-22', 'unc-89', 'myo-2', 'myo-3', 'myo-1', 'unc-15', 'let-805', 'anc-1'],
                'F2: 127N, Control, R1': [103.6, 99.7, 104.9, 100.4, 103.4, 102.3, 89.0, 96.8, 100.3],
                'F2: 127C, Control, R2': [115.1, 116.1, 98.0, 118.6, 107.2, 107.2, 149.9, 109.9, 94.1],
                'F2: 128N, Control, R3': [90.0, 97.8, 94.5, 97.3, 97.3, 88.2, 103.0, 98.5, 100.8],
                'F2: 128C, Control, R4': [90.8, 99.3, 95.4, 98.5, 85.1, 106.7, 101.3, 103.2, 95.0],
                'F2: 129N, Sample, R1': [105.2, 99.0, 104.3, 100.2, 104.3, 93.4, 98.0, 101.1, 102.1],
                'F2: 130N, Sample, R2': [100.0, 96.6, 101.6, 97.5, 103.8, 83.1, 97.6, 100.5, 100.1],
                'F2: 130C, Sample, R3': [95.4, 96.1, 94.1, 97.8, 104.1, 87.4, 99.1, 102.6, 99.5],
                'F2: 131, Sample, R4': [99.9, 102.3, 97.8, 98.1, 105.0, 87.6, 98.7, 97.4, 100.9]
            }
            df = pd.DataFrame(data)
            df.set_index('Gene Symbol', inplace=True)
        else:
            # Determine file type and read accordingly
            file_extension = uploaded_file.name.split('.')[-1].lower()
            
            if file_extension in ['xlsx', 'xls']:
                df = pd.read_excel(uploaded_file, index_col=0)
            elif file_extension == 'csv':
                df = pd.read_csv(uploaded_file, index_col=0)
            elif file_extension in ['tsv', 'txt']:
                df = pd.read_csv(uploaded_file, sep='\t', index_col=0)
            
        # Display raw data
        st.markdown("<h2 class='sub-header'>2. Review your data</h2>", unsafe_allow_html=True)
        st.dataframe(df, use_container_width=True)
        
        # Data preprocessing options
        st.markdown("<h2 class='sub-header'>3. Configure analysis</h2>", unsafe_allow_html=True)
        
        # Handle missing values
        if df.isnull().values.any():
            st.warning("Your data contains missing values. They will be imputed.")
        
        # Extract column metadata
        col_info = {}
        for col in df.columns:
            # Extract sample information
            sample_info = {}
            
            # Try different patterns to extract sample metadata
            patterns = [
                # Pattern like "F2: 127N, Control, R1"
                r"(?:.*?,\s*)?(.*?)(?:,\s*)(R\d+|Rep\d+|Replicate\d+)$",
                # Pattern like "Control_R1" or "Sample_R2"
                r"(.+?)[-_]+(R\d+|Rep\d+|Replicate\d+)$",
                # Check for numerical replicate indicators at the end
                r"(.+?)[-_]+(\d+)$",
                # Try to extract any word followed by any number
                r"(.+?)(\d+)$"
            ]
            
            metadata_extracted = False
            for pattern in patterns:
                match = re.search(pattern, col)
                if match:
                    sample_info['group'] = match.group(1).strip()
                    sample_info['replicate'] = match.group(2).strip()
                    metadata_extracted = True
                    break
            
            # If no pattern matches, use the entire column name
            if not metadata_extracted:
                sample_info['group'] = col
                sample_info['replicate'] = 'NA'
            
            col_info[col] = sample_info
        
        # Identify unique groups
        unique_groups = list(set(info['group'] for info in col_info.values()))
        
        # Group selection
        with st.expander("Sample grouping (automatically detected)", expanded=True):
            st.info("The tool has automatically detected the following sample groups. You can modify if needed.")
            
            # Let users modify sample groups
            group_mapping = {}
            for i, col in enumerate(df.columns):
                detected_group = col_info[col]['group']
                detected_replicate = col_info[col]['replicate']
                
                col1, col2 = st.columns([3, 1])
                with col1:
                    group_mapping[col] = st.selectbox(
                        f"Group for {col}", 
                        unique_groups + ['Custom'],
                        index=unique_groups.index(detected_group) if detected_group in unique_groups else 0,
                        key=f"group_{i}"
                    )
                    
                    if group_mapping[col] == 'Custom':
                        group_mapping[col] = st.text_input(f"Custom group name for {col}", detected_group, key=f"custom_group_{i}")
                
                with col2:
                    st.text_input("Replicate ID", value=detected_replicate, key=f"replicate_{i}", disabled=True)
        
        # Visualization settings
        with st.expander("Visualization settings", expanded=True):
            col1, col2 = st.columns(2)
            with col1:
                pca_dimensions = st.selectbox("PCA dimensions to display", ["2D", "3D"], index=0)
                pca_components = st.slider("Number of PCs to calculate", min_value=2, max_value=min(10, len(df.columns)), value=min(5, len(df.columns)))
                
                # Add ellipse transparency control
                ellipse_transparency = st.slider(
                    "Confidence region transparency", 
                    min_value=0.0, 
                    max_value=1.0, 
                    value=0.15, 
                    step=0.05,
                    help="Control how transparent the confidence ellipses appear (0 = invisible, 1 = solid)"
                )
                
            with col2:
                st.markdown("#### Color Settings")
                color_options = ["Set1", "Dark2", "Paired", "tab10", "viridis", "plasma"]
                color_palette = st.selectbox(
                    "Color palette", 
                    color_options,
                    index=0,
                    help="Select a color scheme for your plots. Changes will apply to all visualizations."
                )
                
                # Visual preview of selected palette
                preview_colors = sns.color_palette(color_palette, n_colors=5)
                preview_colors_hex = [f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}" for r, g, b in preview_colors]
                
                preview_html = '<div style="display: flex; justify-content: space-between; margin-top: 10px; margin-bottom: 15px;">'
                for color in preview_colors_hex:
                    preview_html += f'<div style="background-color: {color}; width: 18%; height: 25px; border-radius: 4px;"></div>'
                preview_html += '</div>'
                st.markdown(preview_html, unsafe_allow_html=True)
                
                marker_size = st.slider("Marker size", min_value=30, max_value=300, value=120)
        
        # Additional options
        with st.expander("Advanced options", expanded=False):
            scaling_method = st.selectbox("Data scaling method", ["Standard (Z-score)", "Min-Max (0-1)", "None"], index=0)
            pc_scaling = st.checkbox("Scale PC coordinates for visualization", value=True)
            
            st.markdown("#### Gene Loading Plot Options")
            include_loading_plot = st.checkbox("Include gene loading plot", value=True)
            
            # Only show these options if loading plot is enabled
            if include_loading_plot:
                col1, col2 = st.columns(2)
                with col1:
                    top_genes_count = st.slider("Number of top genes to highlight", min_value=5, max_value=30, value=10)
                with col2:
                    arrow_thickness = st.slider("Arrow thickness", min_value=1, max_value=5, value=3)
                    
                highlight_style = st.radio(
                    "Gene highlight style",
                    ["Default", "Bold", "Colorful"],
                    horizontal=True,
                    help="Choose how to highlight top contributing genes in the loading plot"
                )
            
        # Run PCA
        if st.button("Run PCA Analysis", type="primary"):
            st.markdown("<h2 class='sub-header'>4. PCA Results</h2>", unsafe_allow_html=True)
            
            # Progress indicator
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            # 1. Prepare the data
            status_text.text("Preparing data...")
            progress_bar.progress(10)
            
            # Get sample metadata
            sample_metadata = pd.DataFrame({
                'Sample': df.columns,
                'Group': [group_mapping.get(col, col_info[col]['group']) for col in df.columns],
                'Replicate': [col_info[col]['replicate'] for col in df.columns]
            })
            
            # Transpose the data for PCA (samples as rows, genes as columns)
            df_pca = df.transpose()
            
            # 2. Handle missing values
            status_text.text("Handling missing values...")
            progress_bar.progress(20)
            
            if df_pca.isnull().values.any():
                imputer = SimpleImputer(strategy='mean')
                df_pca_imputed = pd.DataFrame(
                    imputer.fit_transform(df_pca),
                    index=df_pca.index,
                    columns=df_pca.columns
                )
            else:
                df_pca_imputed = df_pca
            
            # 3. Scale the data
            status_text.text("Scaling data...")
            progress_bar.progress(30)
            
            if scaling_method == "Standard (Z-score)":
                scaler = StandardScaler()
                df_pca_scaled = pd.DataFrame(
                    scaler.fit_transform(df_pca_imputed),
                    index=df_pca_imputed.index,
                    columns=df_pca_imputed.columns
                )
            elif scaling_method == "Min-Max (0-1)":
                df_pca_scaled = (df_pca_imputed - df_pca_imputed.min()) / (df_pca_imputed.max() - df_pca_imputed.min())
            else:
                df_pca_scaled = df_pca_imputed
            
            # 4. Run PCA
            status_text.text("Running PCA...")
            progress_bar.progress(40)
            
            # Calculate PCA
            pca = PCA(n_components=pca_components)
            pca_result = pca.fit_transform(df_pca_scaled)
            
            # Create a DataFrame with PCA results
            pca_df = pd.DataFrame(
                data=pca_result,
                columns=[f'PC{i+1}' for i in range(pca_components)],
                index=df_pca_scaled.index
            )
            
            # Add group information
            pca_df['Group'] = [group_mapping.get(col, col_info[col]['group']) for col in df_pca_scaled.index]
            
            # 5. Create visualizations
            status_text.text("Generating visualizations...")
            progress_bar.progress(50)
            
            # Prepare for plotting
            unique_groups = list(set(pca_df['Group']))
            
            # Force different color palettes based on selection
            if color_palette == 'Set1':
                colors = sns.color_palette('Set1', n_colors=len(unique_groups))
            elif color_palette == 'Dark2':
                colors = sns.color_palette('Dark2', n_colors=len(unique_groups))
            elif color_palette == 'Paired':
                colors = sns.color_palette('Paired', n_colors=len(unique_groups))
            elif color_palette == 'tab10':
                colors = sns.color_palette('tab10', n_colors=len(unique_groups))
            elif color_palette == 'viridis':
                colors = sns.color_palette('viridis', n_colors=len(unique_groups))
            elif color_palette == 'plasma':
                colors = sns.color_palette('plasma', n_colors=len(unique_groups))
            else:
                colors = sns.color_palette('Set1', n_colors=len(unique_groups))
                
            group_color_map = dict(zip(unique_groups, colors))
            
            # Set a different color palette for the variance plot
            if color_palette in ['viridis', 'plasma']:
                variance_colors = sns.color_palette('Blues_d', n_colors=2)
            else:
                variance_colors = ['skyblue', 'red']
            
            # Plot variance explained
            explained_variance = pca.explained_variance_ratio_ * 100
            cumulative_variance = np.cumsum(explained_variance)
            
            progress_bar.progress(60)
            
            # Create figure for variance plot
            fig_var, ax_var = plt.subplots(figsize=(10, 6))
            
            ax_var.bar(range(1, len(explained_variance) + 1), explained_variance, alpha=0.7, 
                    color=variance_colors[0], label='Individual explained variance')
            ax_var.step(range(1, len(cumulative_variance) + 1), cumulative_variance, where='mid', 
                     color=variance_colors[1], label='Cumulative explained variance')
            
            ax_var.set_xlabel('Principal Component')
            ax_var.set_ylabel('Explained Variance (%)')
            ax_var.set_title('Scree Plot - Explained Variance by Principal Components')
            ax_var.set_xticks(range(1, len(explained_variance) + 1))
            ax_var.set_xticklabels([f'PC{i}' for i in range(1, len(explained_variance) + 1)])
            ax_var.grid(True, linestyle='--', alpha=0.7)
            ax_var.legend()
            
            # Show variance plot
            st.subheader("Explained Variance")
            st.pyplot(fig_var)
            
            # Table with variance explained
            variance_df = pd.DataFrame({
                'Principal Component': [f'PC{i+1}' for i in range(len(explained_variance))],
                'Explained Variance (%)': explained_variance,
                'Cumulative Variance (%)': cumulative_variance
            })
            
            st.dataframe(variance_df, use_container_width=True)
            
            progress_bar.progress(70)
            
            # Create PCA plot
            if pca_dimensions == "2D":
                # Create 2D PCA plot
                fig_pca, ax_pca = plt.subplots(figsize=(12, 8))
                
                # Plot each group
                for group in unique_groups:
                    group_data = pca_df[pca_df['Group'] == group]
                    ax_pca.scatter(
                        group_data['PC1'], 
                        group_data['PC2'],
                        s=marker_size,
                        color=group_color_map[group],
                        label=group,
                        alpha=0.7,
                        edgecolors='k'
                    )
                    
                    # Add sample labels
                    for idx in group_data.index:
                        ax_pca.annotate(
                            idx,
                            (group_data.loc[idx, 'PC1'], group_data.loc[idx, 'PC2']),
                            xytext=(5, 5),
                            textcoords='offset points',
                            fontsize=9
                        )
                
                # Plot properties
                ax_pca.set_xlabel(f'PC1 ({explained_variance[0]:.2f}%)')
                ax_pca.set_ylabel(f'PC2 ({explained_variance[1]:.2f}%)')
                ax_pca.set_title('PCA: Principal Component Analysis')
                ax_pca.grid(True, linestyle='--', alpha=0.7)
                ax_pca.legend(title="Sample Groups")
                
                # Add confidence ellipses for each group
                from matplotlib.patches import Ellipse
                import matplotlib.transforms as transforms
                
                def confidence_ellipse(x, y, ax, n_std=2.0, facecolor='none', **kwargs):
                    if len(x) < 3:
                        return None
                        
                    cov = np.cov(x, y)
                    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
                    
                    # Using a special case to obtain the eigenvalues of this
                    # two-dimensional dataset.
                    ell_radius_x = np.sqrt(1 + pearson)
                    ell_radius_y = np.sqrt(1 - pearson)
                    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2, **kwargs)
                    
                    # Calculating the standard deviation of x from the square root of the variance
                    scale_x = np.sqrt(cov[0, 0]) * n_std
                    # Calculating the standard deviation of y from the square root of the variance
                    scale_y = np.sqrt(cov[1, 1]) * n_std
                    
                    # Scaling the ellipse
                    mean_x = np.mean(x)
                    mean_y = np.mean(y)
                    
                    transf = transforms.Affine2D() \
                        .rotate_deg(45) \
                        .scale(scale_x, scale_y) \
                        .translate(mean_x, mean_y)
                    
                    ellipse.set_transform(transf + ax.transData)
                    
                    return ax.add_patch(ellipse)
                
                # Add 95% confidence ellipses - with more transparent fill
                for group in unique_groups:
                    if len(pca_df[pca_df['Group'] == group]) >= 3:  # Need at least 3 points
                        # Get the group color but make it very transparent for the fill
                        group_color_rgb = group_color_map[group]
                        
                        # Draw the ellipse with a very faded fill
                        confidence_ellipse(
                            pca_df[pca_df['Group'] == group]['PC1'],
                            pca_df[pca_df['Group'] == group]['PC2'],
                            ax_pca,
                            n_std=2.0,
                            edgecolor=group_color_map[group],
                            facecolor=group_color_rgb,
                            alpha=ellipse_transparency,  # Use the value from the slider
                            linestyle='--',
                            linewidth=1.5
                        )
                
            else:  # 3D PCA plot
                from mpl_toolkits.mplot3d import Axes3D
                
                fig_pca = plt.figure(figsize=(14, 10))
                ax_pca = fig_pca.add_subplot(111, projection='3d')
                
                # Plot each group
                for group in unique_groups:
                    group_data = pca_df[pca_df['Group'] == group]
                    ax_pca.scatter(
                        group_data['PC1'], 
                        group_data['PC2'],
                        group_data['PC3'],
                        s=marker_size,
                        color=group_color_map[group],
                        label=group,
                        alpha=0.7,
                        edgecolors='k'
                    )
                    
                    # Add sample labels
                    for idx in group_data.index:
                        ax_pca.text(
                            group_data.loc[idx, 'PC1'], 
                            group_data.loc[idx, 'PC2'],
                            group_data.loc[idx, 'PC3'],
                            idx,
                            fontsize=9
                        )
                
                # Plot properties
                ax_pca.set_xlabel(f'PC1 ({explained_variance[0]:.2f}%)')
                ax_pca.set_ylabel(f'PC2 ({explained_variance[1]:.2f}%)')
                ax_pca.set_zlabel(f'PC3 ({explained_variance[2]:.2f}%)')
                ax_pca.set_title('PCA: Principal Component Analysis (3D)')
                ax_pca.legend(title="Sample Groups")
            
            # Show PCA plot
            st.subheader("PCA Plot")
            st.pyplot(fig_pca)
            
            progress_bar.progress(80)
            
            # Create loadings plot if requested
            if include_loading_plot:
                # Get the loadings
                loadings = pca.components_.T
                loadings_df = pd.DataFrame(
                    loadings,
                    columns=[f'PC{i+1}' for i in range(pca_components)],
                    index=df.index
                )
                
                # Calculate the contribution of each gene
                gene_contribution = pd.DataFrame(
                    index=loadings_df.index,
                    columns=['PC1_contribution', 'PC2_contribution', 'Total_contribution']
                )
                
                gene_contribution['PC1_contribution'] = loadings_df['PC1'].abs() * explained_variance[0]
                gene_contribution['PC2_contribution'] = loadings_df['PC2'].abs() * explained_variance[1]
                gene_contribution['Total_contribution'] = gene_contribution['PC1_contribution'] + gene_contribution['PC2_contribution']
                
                # Get top contributing genes
                top_genes = gene_contribution.sort_values('Total_contribution', ascending=False).head(top_genes_count)
                
                # Create loadings plot
                fig_loadings, ax_loadings = plt.subplots(figsize=(12, 12))
                
                # Plot all genes with light color
                ax_loadings.scatter(
                    loadings_df['PC1'],
                    loadings_df['PC2'],
                    s=40,
                    color='lightgray',
                    alpha=0.3
                )
                
                # Determine color schemes based on settings
                if highlight_style == "Colorful":
                    # Use multiple colors from the selected palette for genes
                    gene_colors = sns.color_palette(color_palette, n_colors=min(top_genes_count, 10))
                    # Repeat colors if needed
                    gene_colors = gene_colors * (top_genes_count // 10 + 1)
                    gene_colors = gene_colors[:top_genes_count]
                    
                    # Create a dictionary mapping genes to colors
                    gene_color_dict = dict(zip(top_genes.index, gene_colors))
                    
                    # Common arrow color
                    arrow_color = 'navy' if color_palette in ['Paired', 'Set1', 'tab10'] else 'black'
                else:
                    # Use a single highlight color
                    if color_palette in ['viridis', 'plasma']:
                        highlight_color = sns.color_palette(color_palette, n_colors=1)[0]
                        arrow_color = sns.color_palette('Blues', n_colors=7)[5]
                    elif color_palette == 'Dark2':
                        highlight_color = sns.color_palette('Dark2', n_colors=3)[0]
                        arrow_color = sns.color_palette('Dark2', n_colors=3)[1]
                    elif color_palette == 'Set1':
                        highlight_color = sns.color_palette('Set1', n_colors=3)[0]
                        arrow_color = sns.color_palette('Set1', n_colors=3)[2]
                    else:
                        highlight_color = 'red'
                        arrow_color = 'blue'
                
                # Plot top genes with appropriate colors
                if highlight_style == "Colorful":
                    for i, gene in enumerate(top_genes.index):
                        ax_loadings.scatter(
                            loadings_df.loc[gene, 'PC1'],
                            loadings_df.loc[gene, 'PC2'],
                            s=120,
                            color=gene_color_dict[gene],
                            alpha=0.9,
                            edgecolors='white',
                            linewidths=1.5,
                            zorder=5
                        )
                else:
                    ax_loadings.scatter(
                        loadings_df.loc[top_genes.index, 'PC1'],
                        loadings_df.loc[top_genes.index, 'PC2'],
                        s=120,
                        color=highlight_color,
                        alpha=0.9,
                        edgecolors='white' if highlight_style == "Bold" else 'k',
                        linewidths=1.5 if highlight_style == "Bold" else 1,
                        zorder=5
                    )
                
                # Add gene labels with style based on selection
                for gene in top_genes.index:
                    label_kwargs = {
                        'fontsize': 12 if highlight_style == "Bold" else 11,
                        'fontweight': 'bold',
                    }
                    
                    if highlight_style in ["Bold", "Colorful"]:
                        label_kwargs['bbox'] = dict(
                            boxstyle="round,pad=0.3", 
                            fc="white", 
                            ec="gray" if highlight_style == "Bold" else gene_color_dict.get(gene, "gray") if highlight_style == "Colorful" else "gray",
                            alpha=0.9,
                            linewidth=1.5
                        )
                    
                    ax_loadings.annotate(
                        gene,
                        (loadings_df.loc[gene, 'PC1'], loadings_df.loc[gene, 'PC2']),
                        xytext=(8, 8),
                        textcoords='offset points',
                        **label_kwargs
                    )
                
                # Set arrow thickness based on slider
                actual_arrow_thickness = arrow_thickness * 0.8  # Convert slider value to appropriate thickness
                
                # Add arrows from origin
                for gene in top_genes.index:
                    # Use gene-specific color if using colorful style
                    arrow_fc = gene_color_dict.get(gene, arrow_color) if highlight_style == "Colorful" else arrow_color
                    arrow_ec = gene_color_dict.get(gene, arrow_color) if highlight_style == "Colorful" else arrow_color
                    
                    ax_loadings.arrow(
                        0, 0,
                        loadings_df.loc[gene, 'PC1'],
                        loadings_df.loc[gene, 'PC2'],
                        head_width=0.02 + (0.01 * arrow_thickness / 3),  # Scale head width with thickness
                        head_length=0.03 + (0.01 * arrow_thickness / 3), # Scale head length with thickness
                        fc=arrow_fc,
                        ec=arrow_ec,
                        linewidth=actual_arrow_thickness,
                        alpha=0.8,
                        zorder=4
                    )
                
                # Add circle
                circle = plt.Circle((0, 0), 1, fill=False, linestyle='--', color='gray')
                ax_loadings.add_patch(circle)
                
                # Set axis limits
                max_val = max(loadings_df['PC1'].abs().max(), loadings_df['PC2'].abs().max()) * 1.2
                ax_loadings.set_xlim(-max_val, max_val)
                ax_loadings.set_ylim(-max_val, max_val)
                
                # Add lines at 0
                ax_loadings.axhline(y=0, color='k', linestyle='-', alpha=0.3)
                ax_loadings.axvline(x=0, color='k', linestyle='-', alpha=0.3)
                
                # Plot properties
                ax_loadings.set_xlabel(f'PC1 Loadings ({explained_variance[0]:.2f}%)')
                ax_loadings.set_ylabel(f'PC2 Loadings ({explained_variance[1]:.2f}%)')
                ax_loadings.set_title('PCA Loading Plot: Gene Contributions')
                ax_loadings.grid(True, linestyle='--', alpha=0.7)
                
                # Show loadings plot
                st.subheader("Gene Contributions (Loading Plot)")
                st.pyplot(fig_loadings)
                
                # Table of top genes
                st.markdown("#### Top Contributing Genes")
                gene_contrib_display = top_genes.copy()
                gene_contrib_display['PC1_contribution'] = gene_contrib_display['PC1_contribution'].map(lambda x: f"{x:.2f}%")
                gene_contrib_display['PC2_contribution'] = gene_contrib_display['PC2_contribution'].map(lambda x: f"{x:.2f}%")
                gene_contrib_display['Total_contribution'] = gene_contrib_display['Total_contribution'].map(lambda x: f"{x:.2f}%")
                st.dataframe(gene_contrib_display, use_container_width=True)
            
            progress_bar.progress(90)
            
            # PCA Interpretation
            st.markdown("<h2 class='sub-header'>5. Interpretation</h2>", unsafe_allow_html=True)
            
            # Calculate group separations
            group_centroids = {}
            for group in unique_groups:
                group_data = pca_df[pca_df['Group'] == group]
                group_centroids[group] = {
                    'PC1': group_data['PC1'].mean(),
                    'PC2': group_data['PC2'].mean()
                }
            
            # Calculate distances between group centroids
            import itertools
            group_distances = {}
            for group1, group2 in itertools.combinations(unique_groups, 2):
                distance = np.sqrt(
                    (group_centroids[group1]['PC1'] - group_centroids[group2]['PC1'])**2 +
                    (group_centroids[group1]['PC2'] - group_centroids[group2]['PC2'])**2
                )
                group_distances[(group1, group2)] = distance
            
            # Calculate within-group variation
            group_variations = {}
            for group in unique_groups:
                group_data = pca_df[pca_df['Group'] == group]
                if len(group_data) > 1:  # Need at least 2 samples
                    centroid = (group_centroids[group]['PC1'], group_centroids[group]['PC2'])
                    distances = []
                    for _, row in group_data.iterrows():
                        point = (row['PC1'], row['PC2'])
                        dist = np.sqrt((point[0] - centroid[0])**2 + (point[1] - centroid[1])**2)
                        distances.append(dist)
                    group_variations[group] = np.mean(distances)
                else:
                    group_variations[group] = 0
            
            # Generate interpretation based on PCA results
            with st.markdown("<div class='interpretation'>", unsafe_allow_html=True):
                st.markdown("### PCA Results Interpretation")
                
                # Explained variance interpretation
                st.markdown("#### Explained Variance")
                pc1_var = explained_variance[0]
                pc2_var = explained_variance[1]
                total_var_first2 = pc1_var + pc2_var
                
                if total_var_first2 > 70:
                    st.markdown(f"PC1 and PC2 together explain **{total_var_first2:.2f}%** of the total variance, which is excellent. This indicates that the 2D representation captures most of the important patterns in your data.")
                elif total_var_first2 > 50:
                    st.markdown(f"PC1 and PC2 together explain **{total_var_first2:.2f}%** of the total variance, which is good. This 2D representation captures a substantial portion of the important patterns.")
                else:
                    st.markdown(f"PC1 and PC2 together explain only **{total_var_first2:.2f}%** of the total variance. This suggests that higher dimensions might be important for your data, and you may want to consider examining PC3 and beyond.")
                
                # Group separation interpretation
                st.markdown("#### Group Separation")
                if len(unique_groups) > 1:
                    # Find max and min distances
                    max_dist_pair = max(group_distances.items(), key=lambda x: x[1])
                    min_dist_pair = min(group_distances.items(), key=lambda x: x[1])
                    
                    st.markdown(f"The greatest separation is between **{max_dist_pair[0][0]}** and **{max_dist_pair[0][1]}** (distance = {max_dist_pair[1]:.2f}).")
                    st.markdown(f"The least separation is between **{min_dist_pair[0][0]}** and **{min_dist_pair[0][1]}** (distance = {min_dist_pair[1]:.2f}).")
                    
                    # Interpretation based on distances
                    if max_dist_pair[1] > 3 * min_dist_pair[1]:
                        st.markdown("There is a **strong separation** between some groups in your data, suggesting clear biological differences between these conditions.")
                    elif max_dist_pair[1] > 1.5 * min_dist_pair[1]:
                        st.markdown("There is a **moderate separation** between groups in your data, suggesting some biological differences between conditions.")
                    else:
                        st.markdown("The separation between groups is **relatively small**, suggesting subtle differences between conditions.")
                
                # Within-group variation interpretation
                st.markdown("#### Within-group Variation")
                if len(group_variations) > 0:
                    max_var_group = max(group_variations.items(), key=lambda x: x[1])
                    min_var_group = min(group_variations.items(), key=lambda x: x[1])
                    
                    st.markdown(f"**{max_var_group[0]}** shows the highest within-group variation (average distance = {max_var_group[1]:.2f}).")
                    st.markdown(f"**{min_var_group[0]}** shows the lowest within-group variation (average distance = {min_var_group[1]:.2f}).")
                    
                    # Interpretation of variation
                    high_var_ratio = max_var_group[1] / (min_var_group[1] if min_var_group[1] > 0 else 0.001)
                    if high_var_ratio > 3:
                        st.markdown(f"There is **significantly higher variability** in the {max_var_group[0]} group compared to others, which may indicate heterogeneity within this group or potential outliers.")
                    elif high_var_ratio > 1.5:
                        st.markdown(f"There is **moderately higher variability** in the {max_var_group[0]} group compared to others.")
                    else:
                        st.markdown("All groups show **relatively similar levels of internal variation**, suggesting consistent replication within each condition.")
                
                # Gene contribution interpretation
                if include_loading_plot:
                    st.markdown("#### Gene Contributions")
                    st.markdown(f"The top {top_genes_count} genes contributing to the separation in your data are highlighted in the loading plot above.")
                    
                    # Top gene for PC1
                    top_pc1_gene = loadings_df['PC1'].abs().sort_values(ascending=False).index[0]
                    top_pc1_value = loadings_df.loc[top_pc1_gene, 'PC1']
                    top_pc1_dir = "positive" if top_pc1_value > 0 else "negative"
                    
                    # Top gene for PC2
                    top_pc2_gene = loadings_df['PC2'].abs().sort_values(ascending=False).index[0]
                    top_pc2_value = loadings_df.loc[top_pc2_gene, 'PC2']
                    top_pc2_dir = "positive" if top_pc2_value > 0 else "negative"
                    
                    st.markdown(f"**{top_pc1_gene}** has the strongest influence on PC1 (in the {top_pc1_dir} direction).")
                    st.markdown(f"**{top_pc2_gene}** has the strongest influence on PC2 (in the {top_pc2_dir} direction).")
                    
                    # Check if any genes have strong loadings in both PCs
                    top_both = gene_contribution['Total_contribution'].sort_values(ascending=False).index[0]
                    if top_both != top_pc1_gene and top_both != top_pc2_gene:
                        st.markdown(f"**{top_both}** has strong contributions to both PC1 and PC2, suggesting it may be a key driver of the overall sample separation.")
                
                # Final summary
                st.markdown("#### Summary")
                if total_var_first2 > 60 and len(unique_groups) > 1 and max(group_distances.values()) > 1.0:
                    st.markdown("Your PCA analysis reveals **clear patterns** in your gene expression data, with good separation between sample groups and identifiable gene drivers.")
                elif total_var_first2 > 40:
                    st.markdown("Your PCA analysis shows **moderate patterns** in your gene expression data. There is some separation between groups, but the differences are not extremely strong.")
                else:
                    st.markdown("Your PCA analysis suggests **subtle patterns** in your gene expression data. The differences between groups are not strongly captured in the first two principal components.")
                
                # Recommendations
                st.markdown("#### Recommendations")
                recommendations = []
                
                if total_var_first2 < 50:
                    recommendations.append("Consider examining higher dimensions (PC3+) as they may contain important information.")
                
                if len(unique_groups) > 1 and min(group_distances.values()) < 0.5:
                    recommendations.append("Some groups show minimal separation, suggesting they may be biologically similar under the conditions studied.")
                
                if max_var_group[1] > 2.0:
                    recommendations.append(f"Investigate potential outliers in the {max_var_group[0]} group due to its high internal variation.")
                
                if recommendations:
                    for rec in recommendations:
                        st.markdown(f"- {rec}")
                else:
                    st.markdown("- Your PCA analysis appears robust. The detected patterns are likely to represent genuine biological differences between your sample groups.")
            
            # Download options
            st.markdown("<h2 class='sub-header'>6. Download Results</h2>", unsafe_allow_html=True)
            
            # Create function for downloading plots
            def get_image_download_link(fig, filename, description, format_type):
                # Save figure to a buffer
                buf = io.BytesIO()
                
                if format_type == 'pdf':
                    # Create PDF with all plots
                    with PdfPages(buf) as pdf:
                        # Add PCA plot
                        pdf.savefig(fig_pca)
                        # Add variance plot
                        pdf.savefig(fig_var)
                        # Add loadings plot if available
                        if include_loading_plot:
                            pdf.savefig(fig_loadings)
                else:
                    # Save individual plot
                    dpi = 300 if format_type in ['png', 'jpg'] else 100
                    fig.savefig(buf, format=format_type, dpi=dpi, bbox_inches='tight')
                
                buf.seek(0)
                b64 = base64.b64encode(buf.read()).decode()
                href = f'<a href="data:image/{format_type};base64,{b64}" download="{filename}.{format_type}">{description}</a>'
                return href
            
            # Generate download links
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.markdown(get_image_download_link(fig_pca, "pca_plot", "Download PCA Plot (PNG)", "png"), unsafe_allow_html=True)
                st.markdown(get_image_download_link(fig_var, "variance_plot", "Download Variance Plot (PNG)", "png"), unsafe_allow_html=True)
            
            with col2:
                st.markdown(get_image_download_link(fig_pca, "pca_plot", "Download PCA Plot (JPG)", "jpg"), unsafe_allow_html=True)
                st.markdown(get_image_download_link(fig_var, "variance_plot", "Download Variance Plot (JPG)", "jpg"), unsafe_allow_html=True)
            
            with col3:
                st.markdown(get_image_download_link(fig_pca, "pca_plot", "Download All Plots (PDF)", "pdf"), unsafe_allow_html=True)
            
            # Download data as CSV
            pca_results_csv = pca_df.to_csv().encode('utf-8')
            st.download_button(
                "Download PCA Results as CSV",
                pca_results_csv,
                "pca_results.csv",
                "text/csv",
                key='download-csv'
            )
            
            if include_loading_plot:
                loadings_csv = loadings_df.to_csv().encode('utf-8')
                st.download_button(
                    "Download Gene Loadings as CSV",
                    loadings_csv,
                    "gene_loadings.csv",
                    "text/csv",
                    key='download-loadings-csv'
                )
            
            # Finish progress
            progress_bar.progress(100)
            status_text.empty()
            
            # Show completion message
            st.success("PCA analysis completed successfully!")
            
    except Exception as e:
        st.error(f"An error occurred during analysis: {str(e)}")
else:
    st.info("Please upload a gene expression data file or use the demo data to get started.")

# Footer
st.markdown("---")
st.markdown("""
<div style="text-align: center;">
    <p>Made with ❤️ for Bioinformatics</p>
    <p>Fork this app on GitHub and contribute!</p>
</div>
""", unsafe_allow_html=True)