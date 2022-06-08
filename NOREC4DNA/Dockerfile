FROM python:3.9-slim as builder
#FROM python:3.6-alpine
MAINTAINER Peter Michael Schwarz "peter.schwarz@uni-marburg.de"
# uwsgi-plugin-python3
COPY . /norec4dna
WORKDIR /norec4dna

RUN apt-get update -y \
 #&& apt-get install  nginx build-essential \
 && apt-get install --no-install-recommends -y gcc build-essential ffmpeg \
 && pip3 install -r requirements.txt --no-cache-dir \
 && python setup.py install \
 && apt-get purge -y --auto-remove build-essential \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# squash / reduce size
FROM scratch
COPY --from=builder / /
WORKDIR /norec4dna

#ENTRYPOINT ["python", "demo_raptor_encode.py", "Dorn", "--error_correction=reedsolomon", "--repairsymbols=6", "--asdna"]
ENTRYPOINT ["python", "find_minimum_packets.py"]
#ENTRYPOINT ["python", "optimize_dist_gd.py"]
# OUTPUT TO parallel_RU10/