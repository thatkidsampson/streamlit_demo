#!/bin/bash

PORT=${PORT:-8501}

streamlit run Home.py \
    --server.port=$PORT \
    --server.address=0.0.0.0 \
    --server.enableXsrfProtection=false \
    --server.enableCORS=false
