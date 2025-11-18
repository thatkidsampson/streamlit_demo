FROM python:3.11-slim

#Update pip
RUN pip install --upgrade pip

# Cloud run expected port
ENV PORT=8080
EXPOSE 8080
# Streamlit local development port
EXPOSE 8501

# Disable Streamlit WebSocket compression (fixes Safari issues)
ENV STREAMLIT_SERVER_ENABLE_WEBSOCKET_COMPRESSION=false

WORKDIR /app

COPY . ./
COPY entrypoint.sh /app/

RUN pip install -r requirements.txt
RUN pip install --no-cache-dir -e .
RUN chmod +x /app/entrypoint.sh

ENTRYPOINT ["/app/entrypoint.sh"]