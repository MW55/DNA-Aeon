FROM python:3.8 as builder
MAINTAINER Peter Michael Schwarz "peter.schwarz@uni-marburg.de"

COPY . /DNA_Aeon
WORKDIR /DNA_Aeon

RUN grep -rl '/home/michael/Code/arithmetic_modulator_error_correction' . | xargs sed -i 's/\/home\/michael\/Code\/arithmetic_modulator_error_correction/\/DNA_Aeon/g'

RUN apt-get update -y \
 && apt-get install --no-install-recommends -y python-dev software-properties-common gcc virtualenv build-essential cmake make\
 && python -m pip install numpy \
 && python -m pip install pandas

# build DNA-Aeon
RUN cmake . && make

# setup NOREC4DNA + dependencies
WORKDIR /DNA_Aeon/NOREC4DNA
RUN find /DNA_Aeon/NOREC4DNA -name '*.pyc' -delete
RUN rm -rf /DNA_Aeon/NOREC4DNA/venv && python -m venv venv && . /DNA_Aeon/NOREC4DNA/venv/bin/activate && pip install wheel && pip install -r ../NOREC_requirements.txt && python setup.py install
WORKDIR /DNA_Aeon

RUN apt-get purge -y --auto-remove build-essential \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# squash / reduce size
FROM scratch
COPY --from=builder / /

WORKDIR /DNA_Aeon
#ENTRYPOINT ["python", "compare2.py"]
#ENTRYPOINT ["bash"]
CMD ["bash"]
