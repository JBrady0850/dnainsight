FROM python:3.12-slim

WORKDIR /app

# Install dependencies first (layer caching)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application source
COPY . .

# Build the bundled SNP reference
RUN python data/build_reference.py

# Create runtime directories
RUN mkdir -p uploads reports_output

# Expose Flask port
EXPOSE 5050

# Run the app
CMD ["python", "app.py"]
