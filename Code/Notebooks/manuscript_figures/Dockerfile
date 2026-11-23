# Exact local base: sha256:f3c900ddf3351b437a66b6af21da2927eaed5b0e25265aaae7ebf0a2a92ec027
# This pre-existing lab image is an external environment input; no fonts/assets added.
FROM python-stats:1.0
RUN python -m pip install --no-cache-dir PyMuPDF==1.26.7
ENV PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
WORKDIR /work
ENTRYPOINT ["python"]
