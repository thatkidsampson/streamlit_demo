#!/bin/bash
# Use PORT environment variable if set, otherwise default to 8501 for local development
PORT=${PORT:-8501}

# Set Streamlit server configuration for Cloud Run
export STREAMLIT_SERVER_HEADLESS=true
export STREAMLIT_SERVER_ENABLECORS=true
export STREAMLIT_SERVER_ENABLEXSRFPROTECTION=false

streamlit run Home.py \
  --server.port=$PORT \
  --server.address=0.0.0.0 \
  --server.headless=true \
  --logger.level=info \
  --server.enableCORS=false \
  --server.enableWebsocketCompression=false