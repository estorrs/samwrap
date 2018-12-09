FROM python:3.6-jessie

RUN apt-get update

# get samtools
RUN git clone https://github.com/samtools/htslib
RUN git clone https://github.com/samtools/samtools
RUN (cd /samtools; autoheader; autoconf -Wno-syntax; ./configure; make; make install)

COPY ./requirements.txt /requirements.txt
RUN pip install -r requirements.txt

COPY . /app
WORKDIR /app

CMD /bin/bash
