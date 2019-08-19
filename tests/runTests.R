library(amap)
library(circlize)
library(scales)
library(shiny)

data("sd01_outputXCMS", package="MetCirc")
data("sd02_deconvoluted", package="MetCirc")

BiocGenerics:::testPackage("MetCirc")