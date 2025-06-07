FROM ubuntu:25.10

# Install dependencies
RUN apt-get update && \
    apt-get install -y build-essential g++ make pkg-config libgtkmm-3.0-dev

# Set the working directory
WORKDIR /app

# Copy the project files into the container
COPY . .

# Build the project
RUN make

# Set the default command (optional, runs your binary)
CMD ["./n_body_simulation"]