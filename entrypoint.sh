#!/bin/bash

PORT=${PORT:-8501}

export STREAMLIT_BROWSER_CWS_ORIGIN="*"
export STREAMLIT_SERVER_COOKIE_SECRET=$(openssl rand -hex 16)
export STREAMLIT_GLOBAL_DISABLE_CONFIG_CHECK=true
export STREAMLIT_SERVER_ENABLE_STATIC_SERVE=false

streamlit run Home.py \
    --server.port=$PORT \
    --server.address=0.0.0.0 \
    --server.enableXsrfProtection=false \
    --server.enableCORS=false
