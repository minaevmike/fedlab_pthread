#-------------------------------------------------
#
# Project created by QtCreator 2014-11-16T15:00:49
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = testmike
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++0x -pthread
LIBS += -pthread

SOURCES += main.cpp
