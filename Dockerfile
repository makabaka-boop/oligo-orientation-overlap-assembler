FROM python:3.12-slim

WORKDIR /app
COPY assembler ./assembler

# The service reads one JSON document from stdin (or a file argument) and
# writes the assembly JSON to stdout.
ENTRYPOINT ["python", "-m", "assembler"]
