FROM continuumio/anaconda3:2019.10

MAINTAINER kmayerblackwell kmayerbl@fredhutch.org

RUN apt-get update && apt-get install -y procps && apt-get install -y nano 

COPY bin/flagellin.py
COPY inputs/*

