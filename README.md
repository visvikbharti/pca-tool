# Gene Expression PCA Analysis Tool

A user-friendly, interactive web application for performing Principal Component Analysis (PCA) on gene expression data.

![PCA Tool Demo](https://github.com/visvikbharti/pca-tool/raw/main/demo_images/demo.png)

## 🔬 About This Tool

This application was designed to help biologists analyze gene expression data without requiring programming knowledge. It automatically identifies sample groups and provides comprehensive visualizations and interpretations of PCA results.

## 🚀 Try It Live

You can access the deployed application here: [Gene Expression PCA Tool](https://visvikbharti-pca-tool.streamlit.app/)

## ✨ Features

- **Intuitive Interface**: User-friendly design for biologists without coding experience
- **Flexible Data Import**: Support for CSV, TSV, and Excel files
- **Automatic Sample Detection**: Smart parsing of sample groups and replicates
- **Dynamic Visualization**: Customizable PCA plots with adjustable settings
- **Gene Loading Analysis**: Identify key genes driving differences between samples
- **Automatic Interpretation**: Get insights on variance, group separation, and gene contributions
- **High-Resolution Exports**: Download publication-ready figures in PNG, JPG, and PDF formats

## 📊 Why Use This Tool?

Standard PCA tools often have limitations:
- They may require equal numbers of replicates across sample groups
- They might have rigid requirements for sample names
- Many lack intuitive interpretation of results

This tool addresses these issues with:
1. **Flexible Sample Grouping**: Handles uneven replicates
2. **Smart Group Detection**: Automatically identifies sample groups from column names
3. **Comprehensive Interpretation**: Provides meaningful biological context to your PCA results

## 🔧 Local Installation

### Prerequisites
- Python 3.7+
- pip (Python package installer)

### Setup Instructions

1. Clone this repository:
   ```bash
   git clone https://github.com/visvikbharti/pca-tool.git
   cd pca-tool
   ```

2. Create a virtual environment:
   ```bash
   python -m venv venv
   
   # On Windows
   venv\Scripts\activate
   
   # On macOS/Linux
   source venv/bin/activate
   ```

3. Install the required packages:
   ```bash
   pip install -r requirements.txt
   ```

4. Run the Streamlit app:
   ```bash
   streamlit run pca-streamlit-app.py
   ```

5. Open your web browser and navigate to http://localhost:8501

## 📁 Input Format

Your input file should follow these guidelines:
- First column: Gene identifiers (symbols or IDs)
- Subsequent columns: Expression values for each sample
- Sample names should ideally include group and replicate information (e.g., "Control_R1", "Treatment_R2")
- Supported formats: CSV, TSV, Excel (.xlsx/.xls)

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 📞 Contact

If you have any questions or feedback, please open an issue on this repository.

---

Made with ❤️ for the scientific community