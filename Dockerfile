FROM continuumio/anaconda3:2019.10

MAINTAINER kmayerblackwell kmayerbl@fredhutch.org

RUN apt-get update && apt-get install -y procps && apt-get install -y nano && apt-get -y install gcc && apt-get -y install unzip && apt-get -y install curl && apt-get -y install wget 

COPY bin/flagellin.py
COPY inputs/*

