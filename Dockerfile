FROM ubuntu:latest

RUN apt-get update && \
    apt-get install -y python3 python3-pip libglib2.0-0 && \
    apt-get clean

COPY modeller_10.7-1_arm64.deb /tmp/modeller_10.7-1_arm64.deb

RUN dpkg -i /tmp/modeller_10.7-1_arm64.deb && \
    rm /tmp/modeller_10.7-1_arm64.deb

ARG KEY_MODELLER
RUN sed -i "s/license = 'xxx'/license = '${KEY_MODELLER}'/" /usr/lib/modeller10.7/modlib/modeller/config.py

ENV PYTHONUNBUFFERED=1

COPY requirements.txt /tmp/requirements.txt
RUN pip3 install --no-cache-dir --break-system-packages -r /tmp/requirements.txt && \
    rm /tmp/requirements.txt

COPY src /app/src
COPY data /app/data
COPY main.py /app/main.py

WORKDIR /app

CMD ["python3", "main.py"]

