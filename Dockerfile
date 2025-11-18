FROM python:3.11

EXPOSE 8080
EXPOSE 8501

ENV STREAMLIT_SERVER_ENABLE_WEBSOCKET_COMPRESSION=false, STREAMLIT_SERVER_CACHE_CONTROL=no-store

WORKDIR /app

# Copy entire project including .streamlit/
COPY . ./

COPY entrypoint.sh /app/
RUN chmod +x /app/entrypoint.sh

RUN pip install -r requirements.txt
RUN pip install --no-cache-dir -e .

ENTRYPOINT ["/app/entrypoint.sh"]