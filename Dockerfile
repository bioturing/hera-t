FROM centos:7

RUN yum -y update

RUN yum -y install gcc make cmake git 

RUN git clone https://github.com/bioturing/hera-t
RUN cd hera-t && git checkout origin/adt
RUN cd hera-t && bash build.sh
RUN cp hera-t/hera-T /usr/bin

CMD /hera-t/hera-T
