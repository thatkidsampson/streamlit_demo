FROM Python:3.11

# Cloud run expected port
EXPOSE 8080
# Streamlit local development port
EXPOSE 8501

WORKDIR /app

COPY . ./

RUN pip install -r requirements.txt
RUN pip install --no-cache-dir -e .

CMD streamlit run Home.py --server.port=$PORT --server.address=0.0.0.0