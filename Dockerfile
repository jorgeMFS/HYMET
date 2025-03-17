# Usar uma imagem base com Conda pré-instalado
FROM continuumio/miniconda3:latest

# Instalar dependências do sistema
RUN apt-get update && apt-get install -y \
    perl \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Copiar o arquivo environment.yml para o contêiner
COPY environment.yml /tmp/environment.yml

# Criar o ambiente Conda a partir do arquivo environment.yml
RUN conda env create -f /tmp/environment.yml

# Adicionar o ambiente ao PATH
ENV PATH /opt/conda/envs/hymet_env/bin:$PATH

# Copiar o código da ferramenta para o contêiner
COPY . /workspace

# Definir o diretório de trabalho
WORKDIR /workspace

# Comando padrão ao executar o contêiner
CMD ["bash"] 