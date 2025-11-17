#!/bin/bash
# Use PORT environment variable if set, otherwise default to 8501 for local development
PORT=${PORT:-8501}

# Set Streamlit server configuration for Cloud Run
export STREAMLIT_SERVER_HEADLESS=true
export STREAMLIT_SERVER_ENABLECORS=false
export STREAMLIT_SERVER_ENABLEXSRFPROTECTION=false
export STREAMLIT_SERVER_ENABLEWEBSOCKETCOMPRESSION=false
export STREAMLIT_CLIENT_LOGGER_LEVEL=debug
# Disable caching to prevent Safari issues
export STREAMLIT_CLIENT_SHOW_ERROR_DETAILS=true

#!/bin/bash
PORT=${PORT:-8501}

# Core settings
export STREAMLIT_SERVER_HEADLESS=true
export STREAMLIT_SERVER_PORT=$PORT
export STREAMLIT_SERVER_ADDRESS=0.0.0.0
export STREAMLIT_SERVER_ENABLECORS=false
export STREAMLIT_SERVER_ENABLEXSRFPROTECTION=false
export STREAMLIT_SERVER_ENABLEWEBSOCKETCOMPRESSION=false

# Safari-specific fixes
export STREAMLIT_BROWSER_GATHERUSAGESTATS=false
export STREAMLIT_SERVER_USESERVERTIME=false

streamlit run Home.py