# HYMET:


Here’s the text in English for your GitHub **README.md**:

---

## Installation and Configuration

Follow the steps below to install and configure the **HYMET** (Hybrid Metagenomic Tool).

### 1. Clone the Repository

First, clone the repository to your local environment:

```bash
git clone https://github.com/your-username/your-repository.git
cd your-repository
```

Replace `your-username/your-repository` with the path to your GitHub repository.

---

### 2. Installation with Docker (Recommended)

If you prefer using Docker, follow these steps:

1. **Build the Docker Image**:
   ```bash
   docker build -t hymet .
   ```

2. **Run the Container**:
   ```bash
   docker run -it hymet
   ```

3. **Inside the Container**:
   - The environment will already be set up with all dependencies installed.
   - Run the tool as needed.

---

### 3. Installation with Conda

If you prefer using Conda, follow these steps:

1. **Create the Conda Environment**:
   ```bash
   conda env create -f environment.yml
   ```

2. **Activate the Environment**:
   ```bash
   conda activate hymet_env
   ```

3. **Run the Tool**:
   - You can now run the tool directly in the terminal.

---

### 4. Running the Tool

After installation (with Docker or Conda), you can run the tool with:

```bash
./run_tool.sh
```

---

### 5. Reference Sketched Databases

The databases required to run the tool are available for download on Google Drive:
- **sketch1.msh.gz**: [Download](https://drive.google.com/drive/folders/1YC0N77UUGinFHNbLpbsucu1iXoLAM6lm?usp=share_link)
- **sketch2.msh.gz**: [Download](https://drive.google.com/drive/folders/1YC0N77UUGinFHNbLpbsucu1iXoLAM6lm?usp=share_link)
- **sketch3.msh.gz**: [Download](https://drive.google.com/drive/folders/1YC0N77UUGinFHNbLpbsucu1iXoLAM6lm?usp=share_link)

Download the files and place them in the `data/` directory of the project.

---

### 6. Requirements

- **Operating System**: Linux, macOS, or Windows.
- **Dependencies**:
  - Conda (for local installation)
  - Docker (for containerized usage)
  - Python 3.8 or higher
  - Perl 5 or higher
  - Minimap2
  - Mash

---

### 7. License

This project is licensed under the **GNU GPL v3**. See the [LICENSE](LICENSE) file for details.

---

### 8. Support

For questions or issues, please open an [issue](https://github.com/your-username/your-repository/issues) in the repository.
