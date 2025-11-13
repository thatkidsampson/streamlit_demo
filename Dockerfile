FROM python:3.11

# Cloud run expected port
EXPOSE 8080
# Streamlit local development port
EXPOSE 8501

WORKDIR /app

COPY . ./
COPY entrypoint.sh /app/
RUN chmod +x /app/entrypoint.sh

RUN pip install -r requirements.txt
RUN pip install --no-cache-dir -e .

ENTRYPOINT ["/app/entrypoint.sh"]