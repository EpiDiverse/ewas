  
FROM nfcore/base
LABEL authors="S. Nilay Can" \
      email="nilaycan@bşologie.uni-marburg.de" \
      description="Docker image containing all requirements for the epidiverse/ewas pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/ewas/bin:$PATH