#!/bin/bash
# Use PORT environment variable if set, otherwise default to 8501 for local development
PORT=${PORT:-8501}
streamlit run Home.py --server.port=$PORT --server.address=0.0.0.0
